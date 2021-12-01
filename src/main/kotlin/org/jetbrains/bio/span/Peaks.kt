package org.jetbrains.bio.span

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.ln
import kotlin.math.log10
import kotlin.math.min

/**
 * The peaks are called in three steps.
 *
 * 1) Firstly, an FDR threshold is applied to the posterior state probabilities
 * 2) Consecutive enriched bins are merged into peaks.
 * 3) Finally, each peak is assigned a qvalue score - median qvalue
 *
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 */
fun SpanFitResults.getPeaks(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    cancellableState: CancellableState? = null
): List<Peak> {
    val progress = Progress { title = "Computing peaks" }.bounded(genomeQuery.get().size.toLong())
    fitInfo.prepareScores()
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        val chromosomePeaks = getChromosomePeaks(chromosome, fdr, gap)
        progress.report(1)
        chromosomePeaks
    }
    progress.done()
    return genomeQuery.get().flatMap { map[it] }
}

fun SpanFitResults.getChromosomePeaks(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int
): List<Peak> =
    // Check that we have information for requested chromosome
    if (chromosome.name in fitInfo.chromosomesSizes) {
        getChromosomePeaks(
            logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
            fitInfo.offsets(chromosome),
            chromosome, fdr, gap,
        )
    } else {
        SpanFitResults.LOG.debug(
            "NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}"
        )
        emptyList()
    }

internal fun SpanFitResults.getChromosomePeaks(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
): List<Peak> {
    val qvalues = PEAKS_QVALUES_CACHE.get(this to chromosome) {
        Fdr.qvalidate(logNullMemberships)
    }
    val enrichedBins = BitterSet(logNullMemberships.size).apply {
        0.until(size()).filter { qvalues[it] < fdr }.forEach(::set)
    }
    Peak.LOG.debug("$chromosome: candidate bins ${enrichedBins.cardinality()}/${logNullMemberships.size}")
    val peaks = enrichedBins.aggregate(gap).map { (i, j) ->
        val passedFDR = (i until j).filter { qvalues[it] < fdr }
        val pvalueLogMedian = DoubleArray(passedFDR.size) { logNullMemberships[passedFDR[it]] }.median()
        val qvalueMedian = DoubleArray(passedFDR.size) { qvalues[passedFDR[it]] }.median()
        val start = offsets[i]
        val end = if (j < offsets.size) offsets[j] else chromosome.length
        // Score should be proportional to length of peak and median original q-value
        val score = min(1000.0, (-log10(qvalueMedian) * (1 + ln((end - start).toDouble())))).toInt()
        // Value is either coverage of fold change
        val value = fitInfo.score(ChromosomeRange(start, end, chromosome))
        Peak(
            chromosome, start, end,
            mlogpvalue = -pvalueLogMedian,
            mlogqvalue = -log10(qvalueMedian),
            value = value,
            score = score
        )
    }
    Peak.LOG.debug("$chromosome: peaks ${peaks.size}")
    return peaks
}


/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * So Q-values estimation is superfluous for each parameter combination.
 * Use cache with weak values to avoid memory overflow.
 */
private val PEAKS_QVALUES_CACHE: Cache<Pair<SpanFitResults, Chromosome>, F64Array> =
    CacheBuilder.newBuilder().weakValues().build()


/**
 * Compute the mean of total coverage from i-th to j-th bin for given labels.
 */
internal fun DataFrame.partialMean(from: Int, to: Int, labelArray: List<String> = labels.toList()) =
    labelArray.map { label -> (from until to).sumOf { getAsInt(it, label) }.toDouble() }.average()


/**
 * This method doesn't duplicate array while computing, so array gets *modified*,
 * instead of StatUtils.percentile(doubles, percentile)
 */
internal fun DoubleArray.median(): Double {
    return object : Percentile(50.0) {
        // force Percentile not to copy scores
        override fun getWorkArray(values: DoubleArray?, begin: Int, length: Int) = this@median
    }.evaluate(this)
}