package org.jetbrains.bio.span

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import com.google.common.collect.ComparisonChain
import org.apache.commons.csv.CSVFormat
import org.apache.commons.math3.stat.descriptive.rank.Percentile
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import java.nio.file.Path
import kotlin.math.ln
import kotlin.math.log10
import kotlin.math.min

/**
 * @param value Foldchange or coverage
 * @param mlogpvalue -log10(pvalue)
 * @param mlogqvalue -log10(qvalue)
 * @param score A score between 0 and 1000
 */
data class Peak(val chromosome: Chromosome,
                val startOffset: Int,
                val endOffset: Int,
                val mlogpvalue: Double,
                val mlogqvalue: Double,
                var value: Double = 0.0,
                val score: Int) : Comparable<Peak>, LocationAware {
    val range: Range
        get() = Range(startOffset, endOffset)

    override val location: Location
        get() = Location(startOffset, endOffset, chromosome, Strand.PLUS)

    override fun compareTo(other: Peak) = ComparisonChain.start()
            .compare(chromosome.name, other.chromosome.name)
            .compare(startOffset, other.startOffset)
            .compare(endOffset, other.endOffset)
            .result()


    companion object {
        internal val LOG = LoggerFactory.getLogger(Peak::class.java)
    }
}

/**
 * Compute the mean of total coverage from i-th to j-th bin for given labels.
 */
internal fun DataFrame.partialMean(from: Int, to: Int, labelArray: List<String> = labels.toList()) =
        labelArray.map { label -> (from until to).map { getAsInt(it, label) }.sum().toDouble() }.average()

/**
 * The peaks are called in three steps.
 *
 * 1) Firstly, an FDR threshold is applied to the posterior state probabilities.
 * 2) Consecutive observations corresponding to H_0 rejections are then merged into ranges.
 * 3) Finally, each range is assigned a qvalue score - median qvalue.
 *
 * @param offsets All the data is binarized, so offsets stores positions in base pair.
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 * @param coverageDataFrame is used to compute either coverage or log fold change.
 */
internal fun getChromosomePeaks(
        logNullMemberships: F64Array,
        qvalues: F64Array,
        offsets: IntArray,
        chromosome: Chromosome,
        fdr: Double,
        gap: Int,
        coverageDataFrame: DataFrame? = null): List<Peak> {
    val enrichedBins = BitterSet(logNullMemberships.size)
    0.until(enrichedBins.size()).filter { qvalues[it] < fdr }.forEach(enrichedBins::set)

    val peaks = enrichedBins.aggregate(gap).map { (i, j) ->
        val passedFDR = (i until j).filter { qvalues[it] < fdr }
        val pvalueLogMedian = DoubleArray(passedFDR.size) { logNullMemberships[passedFDR[it]] }.median()
        val qvalueMedian = DoubleArray(passedFDR.size) { qvalues[passedFDR[it]] }.median()
        val start = offsets[i]
        val end = if (j < offsets.size) offsets[j] else chromosome.length
        // Score should be proportional to length of peak and median original q-value
        val score = min(1000.0, (-log10(qvalueMedian) * (1 + ln((end - start).toDouble()))))
        // Value is either coverage of fold change
        var value = 0.0
        if (coverageDataFrame != null) {
            if (coverageDataFrame.labels.size == 1 ||
                coverageDataFrame.labels.all { it.startsWith(SpanPeakCallingExperiment.TRACK_PREFIX) }
            ) {
                value = coverageDataFrame.partialMean(i, j)
            } else if (coverageDataFrame.labels.all {
                    it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK1_PREFIX) ||
                            it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK2_PREFIX)
                }) {
                val track1 = coverageDataFrame.partialMean(i, j, coverageDataFrame.labels
                    .filter {
                        it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK1_PREFIX)
                    })
                val track2 = coverageDataFrame.partialMean(i, j, coverageDataFrame.labels
                    .filter {
                        it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK2_PREFIX)
                    })
                // Value if LogFC
                value = if (track2 != 0.0) ln(track1) - ln(track2) else Double.MAX_VALUE
            } else {
                Peak.LOG.debug("Failed to compute value for ${coverageDataFrame.labels}")
            }
        }
        Peak(
            chromosome, start, end,
            mlogpvalue = -pvalueLogMedian,
            mlogqvalue = -log10(qvalueMedian),
            value = value,
            score = score.toInt()
        )
    }
    return peaks
}


/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * So Q-values estimation is superfluous for each parameters combination.
 * Use cache with weak values to avoid memory overflow.
 */
private val f64QValuesCache: Cache<Pair<SpanFitResults, Chromosome>, F64Array> =
        CacheBuilder.newBuilder().weakValues().build()


fun SpanFitResults.getChromosomePeaks(
        chromosome: Chromosome,
        fdr: Double,
        gap: Int,
        coverageDataFrame: DataFrame? = null
): List<Peak> =
        // Check that we have information for requested chromosome
        if (chromosome.name in fitInfo.chromosomesSizes) {
            val f64LogNullMemberships =
                    logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)
            val f64QValues = f64QValuesCache.get(this to chromosome) {
                Fdr.qvalidate(f64LogNullMemberships)
            }
            getChromosomePeaks(
                    f64LogNullMemberships,
                    f64QValues,
                    fitInfo.offsets(chromosome),
                    chromosome, fdr, gap,
                    coverageDataFrame)
        } else {
            SpanFitResults.LOG.debug(
                    "NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}")
            emptyList()
        }

fun SpanFitResults.getPeaks(
        genomeQuery: GenomeQuery,
        fdr: Double,
        gap: Int,
        cancellableState: CancellableState? = null
): List<Peak> {
    val coverage = fitInfo.scoresDataFrame()
    if (coverage.isEmpty()) {
        SpanFitResults.LOG.debug("No coverage caches present, peak scores won't be computed.")
    }
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        getChromosomePeaks(chromosome, fdr, gap, coverage[chromosome])
    }
    return genomeQuery.get().flatMap { map[it] }
}

fun savePeaks(peaks: List<Peak>, path: Path, peakName: String = "") {
    Peak.LOG.debug("""FORMAT $path:
chromosome, start, end, name, score, strand, coverage/foldchange, -log(pvalue), -log(qvalue)""")
    CSVFormat.TDF.print(path.bufferedWriter()).use { printer ->
        peaks.sorted().forEachIndexed { i, peak ->
            /* See MACS2 output format for details https://github.com/taoliu/MACS/ */
            printer.printRecord(
                    peak.chromosome.name,
                    peak.range.startOffset.toString(),
                    peak.range.endOffset.toString(),
                    "${if (peakName.isNotEmpty()) peakName else "peak"}_${i + 1}",
                    peak.score.toString(),
                    ".",
                    peak.value.toString(),
                    peak.mlogpvalue.toString(),
                    peak.mlogqvalue.toString())
        }
    }
}

/**
 * This method doesn't duplicate array while computing, so array gets *modified*,
 * instead of StatUtils.percentile(doubles, percentile)
 */
fun DoubleArray.median(): Double {
    return object : Percentile(50.0) {
        // force Percentile not to copy scores
        override fun getWorkArray(values: DoubleArray?, begin: Int, length: Int) = this@median
    }.evaluate(this)
}