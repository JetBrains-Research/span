package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.span.coverage.SpanCoverageSampler.sampleCoverage
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_HARD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_SPEED
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_LIGHT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SENSITIVITY
import org.jetbrains.bio.span.peaks.SpanModelToPeaks
import org.jetbrains.bio.util.Logs
import org.jetbrains.bio.util.stemGz
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import org.slf4j.event.Level
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * @author Oleg Shpynov
 * @date 22/11/2018
 */
class SpanModelFitExperimentTest {

    @Test
    fun testDataQuery() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->
            sampleCoverage(path, GenomeQuery(Genome["to1"]), 200, goodQuality = true)
            println("Saved sampled track file: $path")
            val dataQuery = SpanAnalyzeFitInformation.createFitInformation(
                GenomeQuery(Genome["to1"]), listOf(SpanDataPaths(path, null)), null,
                listOf("foo"), AutoFragment, true, 200
            ).dataQuery

            assertTrue(dataQuery.id.startsWith("${path.stemGz}_foo_200"))
            val df = dataQuery.apply(Chromosome(Genome["to1"], "chr1"))
            assertEquals("[foo]", df.labels.toList().toString())
        }
    }


    @Test
    fun testEffectiveGenomeQueryEmptyChromosomes() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->
            sampleCoverage(
                path, GenomeQuery(Genome["to1"], "chr1"), 200, goodQuality = true
            )
            println("Saved sampled track file: $path")
            val (out, _) = Logs.captureLoggingOutput(Level.DEBUG) {
                val effectiveGenomeQuery = SpanModelFitExperiment.filterGenomeQueryWithData(
                    GenomeQuery(Genome["to1"]), listOf(SpanDataPaths(path, null)), null, AutoFragment, true
                )
                assertEquals("[chr1]", effectiveGenomeQuery.get().map { it.name }.toString())
            }
            assertIn("Chromosomes with no reads detected are ignored. Use --debug for details.", out)
            assertIn("Ignored chromosomes: chr2,chr3,chrX,chrM", out)
        }
    }

    @Test
    fun analyzeExtendedChromosomePeaksLoaded() {
        val effectiveGenomeQuery = GenomeQuery(Genome["to1"], "chr1")
        val fullGenomeQuery = GenomeQuery(Genome["to1"])
        val emptyMaps = genomeMap(effectiveGenomeQuery) { BitSet() }
        withTempFile("track", ".bed.gz") { path ->
            sampleCoverage(
                path, effectiveGenomeQuery, 200, emptyMaps, emptyMaps, goodQuality = true
            )
            println("Saved sampled track file: $path")
            val peakCallingExperiment = SpanPeakCallingExperiment.getExperiment(
                fullGenomeQuery, listOf(SpanDataPaths(path, null)), null, AutoFragment, true, 200
            )
            assertTrue(
                SpanModelToPeaks.getPeaks(
                    peakCallingExperiment.results, fullGenomeQuery,
                    fdr = 0.05,
                    multipleTesting = SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION,
                    sensitivityCmdArg = SPAN_DEFAULT_SENSITIVITY,
                    gapCmdArg = SPAN_DEFAULT_GAP,
                    summits = false,
                    fragmentationLight = SPAN_DEFAULT_FRAGMENTATION_LIGHT,
                    fragmentationHard = SPAN_DEFAULT_FRAGMENTATION_HARD,
                    fragmentationSpeed = SPAN_DEFAULT_FRAGMENTATION_SPEED,
                    clip = SPAN_DEFAULT_CLIP_MAX_SIGNAL,
                    cancellableState = null
                )
                    .toList()
                    .isNotEmpty(),
                "Expected peak set not to be empty."
            )
        }
    }

    companion object {
        private fun assertIn(substring: String, fullString: String) {
            // Process Windows with different line separators correctly.
            substring.lines().forEach { s ->
                assertTrue(s in fullString, "Expected <$s> to be in <$fullString>.")
            }
        }
    }
}