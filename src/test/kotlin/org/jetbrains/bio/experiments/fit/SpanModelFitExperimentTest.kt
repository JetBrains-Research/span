package org.jetbrains.bio.experiments.fit

import org.apache.log4j.Level
import org.jetbrains.bio.Logs
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.span.SpanCLALongTest
import org.jetbrains.bio.span.getPeaks
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import java.io.ByteArrayOutputStream
import java.io.PrintStream
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
            SpanCLALongTest.sampleCoverage(path, GenomeQuery("to1"), 200, goodQuality = true)
            println("Saved sampled track file: $path")
            val (_, dataQuery) =
                    SpanModelFitExperiment.createEffectiveQueries(GenomeQuery("to1"),
                            listOf(path to null), listOf("foo"), null, 200)

            assertEquals("${path.stemGz}_200_foo", dataQuery.id)
            val df = dataQuery.apply(Chromosome("to1", "chr1"))
            assertEquals("[foo]", df.labels.toList().toString())
        }
    }


    @Test
    fun testEffectiveGenomeQuery() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->
            SpanCLALongTest.sampleCoverage(path, GenomeQuery("to1"), 200, goodQuality = true)
            println("Saved sampled track file: $path")
            val (effectiveGenomeQuery, _) =
                    SpanModelFitExperiment.createEffectiveQueries(GenomeQuery("to1"),
                            listOf(path to null), listOf("foo"), null, 200)
            assertEquals("[chr1, chr2, chr3, chrX]", effectiveGenomeQuery.get().map { it.name }.toString())
        }
    }

    @Test
    fun testEffectiveGenomeQueryEmptyChromosomes() {
        val out = System.out
        val outStream = ByteArrayOutputStream()
        System.setOut(PrintStream(outStream))
        // Update with changed System.out
        Logs.addConsoleAppender(Level.INFO)

        try {
            // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
            withTempFile("track", ".bed.gz") { path ->
                SpanCLALongTest.sampleCoverage(path, GenomeQuery("to1", "chr1"), 200, goodQuality = true)
                println("Saved sampled track file: $path")
                val (effectiveGenomeQuery, _) =
                        SpanModelFitExperiment.createEffectiveQueries(GenomeQuery("to1"),
                                listOf(path to null), listOf("foo"), null, 200)
                assertEquals("[chr1]", effectiveGenomeQuery.get().map { it.name }.toString())
                val sout = String(outStream.toByteArray())
                assertTrue("chr2: no reads detected, ignoring" in sout)
                assertTrue("chr3: no reads detected, ignoring" in sout)
                assertTrue("chrX: no reads detected, ignoring" in sout)
            }
        } finally {
            System.setOut(out)
        }
    }

    @Test
    fun analyzeExtendedChromosomePeaksLoaded() {
        val effectiveGenomeQuery = GenomeQuery("to1", "chr1")
        val fullGenomeQuery = GenomeQuery("to1")
        val emptyMaps = genomeMap(effectiveGenomeQuery) { BitSet() }
        withTempFile("track", ".bed.gz") { path ->
            SpanCLALongTest.sampleCoverage(path, effectiveGenomeQuery, 200, emptyMaps, emptyMaps, goodQuality = true)
            println("Saved sampled track file: $path")
            val peakCallingExperiment =
                    SpanPeakCallingExperiment.getExperiment(fullGenomeQuery, listOf(path to null), 200, null)
            assertTrue(peakCallingExperiment.results.getPeaks(fullGenomeQuery, 0.05, 0).isNotEmpty())
        }
    }
}