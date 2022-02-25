package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.span.SpanCLALongTest
import org.jetbrains.bio.span.peaks.getPeaks
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
            SpanCLALongTest.sampleCoverage(path, GenomeQuery(Genome["to1"]), 200, goodQuality = true)
            println("Saved sampled track file: $path")
            val dataQuery = SpanAnalyzeFitInformation.createFitInformation(
                GenomeQuery(Genome["to1"]), listOf(SpanDataPaths(path, null)),
                listOf("foo"), AutoFragment, true, 200
            ).dataQuery

            assertTrue(dataQuery.id.startsWith("binned_200_${path.stemGz}_foo#"))
            val df = dataQuery.apply(Chromosome(Genome["to1"], "chr1"))
            assertEquals("[foo]", df.labels.toList().toString())
        }
    }


    @Test
    fun testEffectiveGenomeQueryEmptyChromosomes() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->
            SpanCLALongTest.sampleCoverage(
                path, GenomeQuery(Genome["to1"], "chr1"), 200, goodQuality = true
            )
            println("Saved sampled track file: $path")
            val (out, _) = Logs.captureLoggingOutput(Level.DEBUG) {
                val effectiveGenomeQuery = SpanModelFitExperiment.filterGenomeQueryWithData(
                    GenomeQuery(Genome["to1"]), listOf(SpanDataPaths(path, null)), AutoFragment, true
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
            SpanCLALongTest.sampleCoverage(
                path, effectiveGenomeQuery, 200, emptyMaps, emptyMaps, goodQuality = true
            )
            println("Saved sampled track file: $path")
            val peakCallingExperiment = SpanPeakCallingExperiment.getExperiment(
                fullGenomeQuery, listOf(SpanDataPaths(path, null)), 200, AutoFragment
            )
            assertTrue(
                peakCallingExperiment.results.getPeaks(fullGenomeQuery, 0.05, 0).isNotEmpty(),
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