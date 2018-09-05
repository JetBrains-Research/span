package org.jetbrains.bio.span

import com.google.common.collect.ComparisonChain
import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Logger
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.Fdr
import org.jetbrains.bio.statistics.data.BitterSet
import org.jetbrains.bio.statistics.data.DataFrame
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.viktor.F64Array
import java.nio.file.Path

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
        internal val LOG = Logger.getLogger(Peak::class.java)
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
 * 1. Firstly, an FDR threshold [fdr] is applied to the posterior state
 *    probabilities.
 * 2. Consecutive observations corresponding to H_0 rejections are then
 *    merged into ranges using [gap].
 * 3. Finally, each range is assigned a qvalue score - median qvalue.
 *    Silly? Maybe.
 *
 * [coverageDataFrame] is used to compute either coverage or log fold change.
 */
internal fun getChromosomePeaks(logNullMemberships: F64Array,
                                offsets: IntArray,
                                chromosome: Chromosome,
                                fdr: Double,
                                gap: Int,
                                coverageDataFrame: DataFrame? = null): List<Peak> {
    // Filter by qvalues
    val qvalues = Fdr.qvalidate(logNullMemberships)
    val enrichedBins = BitterSet(logNullMemberships.size)
    0.until(enrichedBins.size()).filter { qvalues[it] < fdr }.forEach(enrichedBins::set)

    return enrichedBins.aggregate(gap).map { (i, j) ->
        val set = (i until j).filter { qvalues[it] < fdr }
        val pvalue = set.map { logNullMemberships[it] }.median()
        val qvalue = set.map { qvalues[it] }.median()
        val start = offsets[i]
        val end = if (j < offsets.size) offsets[j] else chromosome.length
        // Score should be proportional to length of peak and average original q-value
        val score = Math.min(1000.0, (-Math.log10(qvalue) * (1 + Math.log((end - start).toDouble()))))
        // Value is either coverage of fold change
        var value = 0.0
        if (coverageDataFrame != null) {
            if (coverageDataFrame.labels.size == 1 ||
                    coverageDataFrame.labels.all { it.startsWith(SpanPeakCallingExperiment.TRACK_PREFIX) }) {
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
                            it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK1_PREFIX)
                        })
                // Value if LogFC
                value = if (track2 != 0.0) Math.log(track1) - Math.log(track2) else Double.MAX_VALUE
            } else {
                Peak.LOG.debug("Failed to compute value for ${coverageDataFrame.labels}")
            }
        }
        Peak(chromosome, start, end,
                mlogpvalue = -pvalue,
                mlogqvalue = -Math.log10(qvalue),
                value = value,
                score = score.toInt())
    }
}

fun SpanFitResults.getChromosomePeaks(chromosome: Chromosome,
                                      fdr: Double,
                                      gap: Int,
                                      coverageDataFrame: DataFrame? = null): List<Peak> =
        getChromosomePeaks(logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
                fitInfo.offsets(chromosome),
                chromosome, fdr, gap,
                coverageDataFrame)

fun SpanFitResults.getPeaks(genomeQuery: GenomeQuery,
                            fdr: Double,
                            gap: Int): List<Peak> {
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        getChromosomePeaks(chromosome, fdr, gap)
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

fun List<Double>.median(): Double {
    if (isEmpty()) throw IllegalStateException("Can't compute median for an empty list.")
    return sorted().let {
        if (it.size % 2 != 0)
            it[it.size / 2]
        else
            (it[it.size / 2 - 1] + it[it.size / 2]) * 0.5
    }
}