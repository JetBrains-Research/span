package org.jetbrains.bio.experiments.histones

import org.apache.commons.math3.stat.StatUtils
import org.apache.log4j.Logger
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.coverage.CoverageWithControl
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.sampling.XSRandom
import org.jetbrains.bio.genome.sequence.BinaryLut
import org.jetbrains.bio.query.BinnedReadsQuery
import org.jetbrains.bio.query.InputQuery
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.util.isAccessible
import org.jetbrains.bio.util.presentablePath
import org.jetbrains.bio.util.size
import java.net.URI
import java.util.stream.Collectors
import java.util.stream.Stream

/**
 * @author Oleg Shpynov
 * @since 11/04/2018.
 */

object PeaksInfo {

    private val LOG = Logger.getLogger(PeaksInfo::class.java)

    private fun Long.formatLongNumber() = String.format("%,d", this).replace(',', ' ')

    fun aboutText(genomeQuery: GenomeQuery,
                  peaksStream: Stream<Location>,
                  src: URI?,
                  readQueries: List<InputQuery<Coverage>>): String {
        val peaks = peaksStream.collect(Collectors.toList())
        val peaksLengths = peaks.map { it.length().toDouble() }.toDoubleArray()
        val peaksCount = peaksLengths.count()
        val peaksLenSum = peaksLengths.sum()
        val coverage = peaksLenSum / genomeQuery.get().map { it.length.toLong() }.sum()
        val srcBlock = if (src != null) """
                    |Track source: ${src.presentablePath()}
                    |Source size: ${src.size}${if (!src.isAccessible()) " <not accessible>" else ""}
                    |
                """.trimMargin() else ""
        val peaksLengthsBlock = """
                    |Peaks Statistics:
                    |  Peaks number: ${peaksCount.toLong().formatLongNumber()}
                    |  Peaks summary length : ${peaksLenSum.toLong().formatLongNumber()} bp
                    |  Genome (${genomeQuery.description}) coverage: ${String.format("%.2f%%", coverage)}
                    |  Lengths:
                    |    Min: ${(if (peaksLengths.isEmpty()) 0L else StatUtils.min(peaksLengths).toLong()).formatLongNumber()} bp
                    |    Mean: ${(if (peaksLengths.isEmpty()) 0L else StatUtils.mean(peaksLengths).toLong()).formatLongNumber()} bp
                    |    Max: ${(if (peaksLengths.isEmpty()) 0L else StatUtils.max(peaksLengths).toLong()).formatLongNumber()} bp
                    |    5%: ${(if (peaksLengths.isEmpty()) 0L else StatUtils.percentile(peaksLengths, 5.0).toLong()).formatLongNumber()} bp
                    |    50%: ${(if (peaksLengths.isEmpty()) 0L else StatUtils.percentile(peaksLengths, 50.0).toLong()).formatLongNumber()} bp
                    |    95%: ${(if (peaksLengths.isEmpty()) 0L else StatUtils.percentile(peaksLengths, 95.0).toLong()).formatLongNumber()} bp
                    |""".trimMargin()
        var signalBlock = ""
        // Don't recompute coverage if it is not processed locally
        if (readQueries.isNotEmpty() && readQueries.all(::cacheAccessible)) {
            val coverages = readQueries.map { it.get() }
            val frip = frip(genomeQuery, peaks, coverages)
            signalBlock += "FRIP: $frip\n"
            val signalToNoise = signalToNoise(genomeQuery, peaks, coverages)
            if (signalToNoise != null) {
                signalBlock += "Signal to noise: $signalToNoise\n"
            }
        }
        return srcBlock + peaksLengthsBlock + signalBlock
    }

    private fun cacheAccessible(query: InputQuery<Coverage>) =
            (query is ReadsQuery && query.npzPath().isAccessible()) ||
                    (query is BinnedReadsQuery && query.bwPath().isAccessible())


    private fun frip(genomeQuery: GenomeQuery, peakLocations: List<Location>, coverages: List<Coverage>): Double {
        val frip = coverages.map(Coverage::signalCoverage).map { coverage ->
            1.0 * peakLocations.map { coverage.getBothStrandsCoverage(it.toChromosomeRange()).toLong() }.sum() /
                    genomeQuery.get().map {
                        coverage.getBothStrandsCoverage(ChromosomeRange(0, it.length, it)).toLong()
                    }.sum()
        }.average()
        LOG.debug("Frip: $frip")
        return frip
    }


    private fun signalToNoise(genomeQuery: GenomeQuery,
                              peaks: List<Location>,
                              coverages: List<Coverage>): Double? {
        val signalCoverages = coverages.map { it.signalCoverage }
        LOG.debug("Generating ${NUMBER_OF_RANGES}x${RANGE_LENGTH}bp top summits")
        val topSummits = topSummits(peaks.toList(), signalCoverages)
        val allPeaksCoords = genomeMap(genomeQuery) { chr ->
            val coords = peaks.filter { it.chromosome == chr }
                    .flatMap { listOf(it.startOffset, it.endOffset) }
                    .sorted()
                    .distinct().toIntArray()
            coords to BinaryLut.of(coords, 24)
        }

        LOG.debug("Generating ${NUMBER_OF_RANGES}x${RANGE_LENGTH}bp desert regions " +
                "with distance $DESERT_DISTANCE")
        val desertRanges = topSummits.mapNotNull {
            sampleRegionInDesert(it.chromosome, allPeaksCoords[it.chromosome])
        }
        if (desertRanges.isEmpty()) {
            LOG.error("Failed to estimate signal-to-noise ratio.")
            return null
        }
        val topPeaksSummitsCoverage = topSummits.map {
            signalCoverages.map { coverage -> coverage.getBothStrandsCoverage(it) }.average()
        }.average()
        val desertCoverage = desertRanges.map {
            signalCoverages.map { coverage -> coverage.getBothStrandsCoverage(it) }.average()
        }.average()
        val signalToNoise = 1.0 * (topPeaksSummitsCoverage + 1e-6) / (desertCoverage + 1e-6)
        LOG.debug("Signal-to-noise ratio $signalToNoise")
        return signalToNoise
    }


    private const val NUMBER_OF_RANGES = 1000
    private const val RANGE_LENGTH = 1000
    private const val DESERT_DISTANCE = 1000

    private val RANDOM = XSRandom()
    private const val SAMPLE_ATTEMPTS_THRESHOLD = 1000

    /**
     * Generates list of given size containing random locations (which may intersect)
     * of given length within chromosome[leftBound, rightBound)
     */
    private fun sampleRegionInDesert(chromosome: Chromosome, allPeaksCoords: Pair<IntArray, BinaryLut>): ChromosomeRange? {
        val (coords, lut) = allPeaksCoords
        var range: ChromosomeRange
        var attempts = 0
        while (true) {
            if (attempts++ > SAMPLE_ATTEMPTS_THRESHOLD) {
                return null
            }
            val startOffset = RANDOM.nextInt(chromosome.length - RANGE_LENGTH)
            range = ChromosomeRange(startOffset, startOffset + RANGE_LENGTH, chromosome)
            if (checkDistance(lut, coords, range.startOffset) && checkDistance(lut, coords, range.endOffset)) {
                return range
            }
        }
    }

    // Check if we are not too close to real peaks
    private fun checkDistance(lut: BinaryLut, coords: IntArray, boundary: Int): Boolean {
        val (low, high) = lut.nearestElemDist(coords, boundary)
        check(low != -1 && high != -1)
        return Math.abs(boundary - coords[low]) > DESERT_DISTANCE && Math.abs(boundary - coords[high]) > DESERT_DISTANCE
    }

    private fun topSummits(peaks: List<Location>, coverages: List<Coverage>): List<ChromosomeRange> {
        return peaks
                .sortedByDescending {
                    coverages.map { coverage ->
                        coverage.getBothStrandsCoverage(it.toChromosomeRange())
                    }.average()
                }
                .map {
                    val start = (it.startOffset + it.endOffset) / 2 - RANGE_LENGTH / 2
                    val end = (it.startOffset + it.endOffset) / 2 + RANGE_LENGTH / 2
                    ChromosomeRange(start, end, it.chromosome)
                }
                .take(NUMBER_OF_RANGES)
                .toList()
    }

}

val Coverage.signalCoverage: Coverage get() = if (this is CoverageWithControl) this.conditionCoverage else this
