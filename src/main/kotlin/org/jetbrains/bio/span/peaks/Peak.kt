package org.jetbrains.bio.span.peaks

import com.google.common.collect.ComparisonChain
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.util.bufferedWriter
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * @param value Foldchange or coverage
 * @param mlogpvalue -log10(pvalue)
 * @param mlogqvalue -log10(qvalue)
 * @param score A score between 0 and 1000
 */
data class Peak(
    val chromosome: Chromosome,
    val startOffset: Int,
    val endOffset: Int,
    val mlogpvalue: Double,
    val mlogqvalue: Double,
    var value: Double = 0.0,
    val score: Int
) : Comparable<Peak>, LocationAware {
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

        fun savePeaks(peaks: List<Peak>, path: Path, peakName: String = "") {
            LOG.debug(
                """FORMAT $path:
chromosome, start, end, name, score, strand, coverage/foldchange, -log(pvalue), -log(qvalue)"""
            )
            CSVFormat.TDF.print(path.bufferedWriter()).use { printer ->
                peaks.sorted().forEachIndexed { i, peak ->
                    /* See MACS2 output format for details https://github.com/taoliu/MACS/ */
                    printer.printRecord(
                        peak.chromosome.name,
                        peak.range.startOffset.toString(),
                        peak.range.endOffset.toString(),
                        "${peakName.ifEmpty { "peak" }}_${i + 1}",
                        peak.score.toString(),
                        ".",
                        peak.value.toString(),
                        peak.mlogpvalue.toString(),
                        peak.mlogqvalue.toString()
                    )
                }
            }
        }
    }
}