package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.span.SpanCLALongTest
import org.jetbrains.bio.span.getPeaks
import org.jetbrains.bio.util.Logs
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
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
            val (_, dataQuery) = SpanModelFitExperiment.createEffectiveQueries(
                GenomeQuery(Genome["to1"]), listOf(path to null),
                listOf("foo"), Optional.empty(), 200
            )

            assertTrue(dataQuery.id.startsWith("${path.stemGz}_200_foo#"))
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
            val (out, _) = Logs.captureLoggingOutput {
                val (effectiveGenomeQuery, _) =
                        SpanModelFitExperiment.createEffectiveQueries(
                            GenomeQuery(Genome["to1"]),
                            listOf(path to null), listOf("foo"), Optional.empty(), 200
                        )
                assertEquals("[chr1]", effectiveGenomeQuery.get().map { it.name }.toString())
            }
            assertIn("chr2: no reads detected, ignoring", out)
            assertIn("chr3: no reads detected, ignoring", out)
            assertIn("chrX: no reads detected, ignoring", out)
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
                fullGenomeQuery, listOf(path to null), 200, Optional.empty()
            )
            assertTrue(
                peakCallingExperiment.results.getPeaks(fullGenomeQuery, 0.05, 0).isNotEmpty(),
                "Expected peak set not to be empty.")
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