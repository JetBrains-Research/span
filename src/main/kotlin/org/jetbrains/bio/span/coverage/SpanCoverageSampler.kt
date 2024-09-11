package org.jetbrains.bio.span.coverage

import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.util.withResource
import java.nio.file.Path
import java.util.*

object SpanCoverageSampler {
    /**
     * Samples coverage from SPAN model on bins resolution
     * @param path coverage file path
     * @param genomeQuery Genome query to sample coverage on
     * @param high bins which should have the highest coverage
     * @param zero bins which should have zero coverage
     * @param goodQuality true if sampling is done from good quality data model
     * @param totalCoverage summary number of reads to sample, null to ignore
     * @param firstZeros Artificially add zeroes to avoid problem with enriched first bin
     */
    fun sampleCoverage(
        path: Path,
        genomeQuery: GenomeQuery, bin: Int,
        high: GenomeMap<BitSet> = genomeMap(genomeQuery) { BitSet() },
        zero: GenomeMap<BitSet> = genomeMap(genomeQuery) { BitSet() },
        goodQuality: Boolean = true,
        totalCoverage: Long? = null,
        firstZeros: Int = 10
    ) {
        withResource(
            SpanCoverageSampler::class.java,
            if (goodQuality)
                "GSM646345_H1_H3K4me3_rep1_hg19_model.json"
            else
                "yd6_k27ac_failed_model.json"
        ) { modelPath ->
            val model = ClassificationModel.load<NB2ZHMM>(modelPath)
            model.logPriorProbabilities[0] = Double.NEGATIVE_INFINITY
            val chromosomes = genomeQuery.get()
            BedFormat().print(path).use { printer ->
                val genomeLength = chromosomes.sumOf { it.length.toLong() }
                chromosomes.forEach { chr ->
                    val bins = chr.length / bin
                    val sampledCoverage = run {
                        repeat(10) {
                            val sample = model.sample(bins).sliceAsInt("d0")
                            if (sample.any { it != 0 }) return@run sample
                        }
                        throw IllegalStateException("Couldn't sample coverage.")
                    }

                    val totalSampledCoverage: Long = sampledCoverage.fold(0L) { acc: Long, v: Int -> acc + v }
                    val scale = if (totalCoverage != null)
                        (chr.length.toFloat() / genomeLength) * (totalCoverage.toFloat() / totalSampledCoverage)
                    else
                        null
                    for (b in firstZeros until bins) {
                        val binIsEnriched = high[chr][b]
                        val binIsZero = zero[chr][b]
                        check(!(binIsEnriched && binIsZero)) { "Both enriched and zero for $b" }
                        val readsInBin = when {
                            binIsEnriched -> bin
                            binIsZero -> 0
                            scale != null -> (scale * sampledCoverage[b]).toInt()
                            else -> sampledCoverage[b]
                        }
                        // Keep reads sorted
                        (0 until readsInBin).map { it % bin }.sorted().forEach { i ->
                            val start = b * bin + i
                            printer.print(BedEntry(chr.name, start, start + 1))
                        }
                    }
                }
            }
        }
    }
}