package org.jetbrains.bio.span.peaks

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanFitResults
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment
import org.jetbrains.bio.experiments.fit.f64Array
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.ln
import kotlin.math.min

const val PEAKS_TYPE_FDR_GAP = "simple"

/**
 * The islands are called in three steps.
 * 1) First, FDR threshold is applied to bins posterior error probabilities.
 * 2) Consecutive enriched bins are merged into peaks.
 * 3) Finally, each peak is assigned a pvalue and qvalue score - median score * ln (length).
 */
fun SpanFitResults.getFdrGapPeaks(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    cancellableState: CancellableState? = null
): List<Peak> {
    val progress = Progress { title = "Computing peaks" }.bounded(genomeQuery.get().size.toLong())
    fitInfo.prepareScores()
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        val chromosomePeaks = getChromosomeFdrGapPeaks(chromosome, fdr, gap)
        progress.report(1)
        chromosomePeaks
    }
    progress.done()
    return genomeQuery.get().flatMap { map[it] }
}

fun SpanFitResults.getChromosomeFdrGapPeaks(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int
): List<Peak> =
    // Check that we have information for requested chromosome
    if (chromosome.name in fitInfo.chromosomesSizes) {
        getChromosomeFdrGapPeaks(
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

internal fun SpanFitResults.getChromosomeFdrGapPeaks(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
): List<Peak> {
    val logFdr = ln(fdr)
    val logQValues = PEAKS_LOG_QVALUES_CACHE.get(this to chromosome) {
        Fdr.qvalidate(logNullMemberships, logResults = true)
    }
    val enrichedBins = BitterSet(logNullMemberships.size).apply {
        0.until(size()).filter { logQValues[it] <= logFdr }.forEach(::set)
    }
    Peak.LOG.debug("$chromosome: candidate bins ${enrichedBins.cardinality()}/${logNullMemberships.size}")
    val peaks = enrichedBins.aggregate(gap).map { (i, j) ->
        val passedFdrBins = (i until j).filter { logQValues[it] <= logFdr }
        val pvalueLogMedian = DoubleArray(passedFdrBins.size) { logNullMemberships[passedFdrBins[it]] }.median()
        val qvalueLogMedian = DoubleArray(passedFdrBins.size) { logQValues[passedFdrBins[it]] }.median()
        val start = offsets[i]
        val end = if (j < offsets.size) offsets[j] else chromosome.length
        // Score should be proportional to length of peak and median original q-value
        val score = min(1000.0, (-qvalueLogMedian * (1 + ln((end - start).toDouble())))).toInt()
        // Value is either coverage of fold change
        val value = fitInfo.score(ChromosomeRange(start, end, chromosome))
        Peak(
            chromosome, start, end,
            mlogpvalue = -pvalueLogMedian / LOG_10,
            mlogqvalue = -qvalueLogMedian / LOG_10,
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
private val PEAKS_LOG_QVALUES_CACHE: Cache<Pair<SpanFitResults, Chromosome>, F64Array> =
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

val LOG_10 = ln(10.0)