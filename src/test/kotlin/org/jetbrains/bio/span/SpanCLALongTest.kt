package org.jetbrains.bio.span

import kotlinx.support.jdk7.use
import org.apache.commons.csv.CSVFormat
import org.apache.log4j.Level
import org.apache.log4j.LogManager
import org.apache.log4j.Logger
import org.jetbrains.bio.Logs
import org.jetbrains.bio.experiments.fit.sampleCoverage
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.ByteArrayOutputStream
import java.io.PrintStream
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class SpanCLALongTest {


    @Before
    fun setUp() {
        (LogManager.getCurrentLoggers().toList() + listOf(LogManager.getRootLogger())).forEach {
            (it as Logger).level = Level.INFO
        }
    }

    @After
    fun tearDown() {
        if (LogManager.getRootLogger().getAppender(Logs.CONSOLE_APPENDER) != null) {
            LogManager.getRootLogger().removeAppender(Logs.CONSOLE_APPENDER)
        }
        System.setOut(OUT)
        System.setErr(ERR)
    }

    @Test
    fun emptyArgs() {
        val stream = ByteArrayOutputStream()
        System.setErr(PrintStream(stream))
        SpanCLA.main(arrayOf())
        assertLinesEqual("""
ERROR: No command given; analyze or compare expected.

Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode (experimental, use at your own risk)
""".trim(), String(stream.toByteArray()).trim())
    }

    @Test
    fun illegalArgs() {
        val stream = ByteArrayOutputStream()
        System.setErr(PrintStream(stream))
        SpanCLA.main(arrayOf("foobar"))
        assertLinesEqual("""
ERROR: Unknown command: foobar; analyze or compare expected.

Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode (experimental, use at your own risk)
""".trim(), String(stream.toByteArray()).trim())
    }

    @Test
    fun quietError() {
        val stream = ByteArrayOutputStream()
        System.setErr(PrintStream(stream))
        SpanCLA.main(arrayOf("foobar", "quiet"))
        assertLinesEqual("""
ERROR: Unknown command: foobar; analyze or compare expected.

Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode (experimental, use at your own risk)
        """.trim(), String(stream.toByteArray()).trim())
    }

    @Test
    fun checkHelp() {
        val stream = ByteArrayOutputStream()
        System.setOut(PrintStream(stream))
        SpanCLA.main(arrayOf("--help"))
        assertLinesEqual("""
Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode (experimental, use at your own risk)
""".trim(), String(stream.toByteArray()).trim())
    }


    @Test
    fun checkVersion() {
        val stream = ByteArrayOutputStream()
        System.setOut(PrintStream(stream))
        SpanCLA.main(arrayOf("--version"))
        assertEquals("@VERSION@.@BUILD@ built on @DATE@", String(stream.toByteArray()).trim())
    }


    // sample random coverage and test the same model prediction at least.
    @Test
    fun compareSameTestOrganismTracks() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val bedPath = it / "peaks.bed"
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                SpanCLA.main(arrayOf("compare",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t1", path.toString(),
                        "-t2", path.toString(),
                        "--output", bedPath.toString(),
                        "--fdr", FDR.toString(),
                        "--gap", GAP.toString(),
                        "--threads", THREADS.toString()))

                FORMAT.parse(bedPath.bufferedReader()).use {
                    for (record in it) {
                        assertTrue(record[7].toDouble() >= Math.log(0.5))
                    }
                }

                val out = String(stream.toByteArray())
                assertIn("""SPAN
COMMAND:
LOG:
WORKING DIR: $it
THREADS: $THREADS
TREATMENT1: $path
CONTROL1: none
TREATMENT2: $path
CONTROL2: none
CHROM.SIZES: $chromsizes
GENOME: to1
FRAGMENT: null
BIN: $BIN
FDR: $FDR
GAP: $GAP
OUTPUT: $bedPath
""", out)
                assertIn("Saved result to $bedPath", out)
                // Check model fit has a progress
                assertIn("] 0.00% (0/250), Elapsed time", out)
            }
        }
    }

    // sample random coverage and test the same model prediction at least.
    @Test
    fun compareSameTestOrganismTracksReplicates() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val bedPath = it / "peaks.bed"
                SpanCLA.main(arrayOf("compare",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", it.toString(),
                        "-b", BIN.toString(),
                        "-g", GAP.toString(),
                        "-fragment", FRAGMENT.toString(),
                        "-t1", "$path,$path",
                        "-t2", "$path,$path,$path",
                        "-o", bedPath.toString(),
                        "--fdr", FDR.toString()))
                FORMAT.parse(bedPath.bufferedReader()).use {
                    for (record in it) {
                        assertTrue(record[7].toDouble() >= 0.5)
                    }
                }
            }
        }
    }

    // sample random coverage and test the same model prediction at least.
    @Test
    fun testModelFitTar() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                SpanCLA.main(arrayOf("analyze",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString()))
                val out = String(stream.toByteArray())
                assertIn("""NO output path given, process model fitting only.
LABELS, FDR, GAP options are ignored.
""", out)
                assertIn(".tar: done in ", out)
                assertIn("Model saved: ", out)
                assertFalse("Loading model" in out)
                assertFalse("Completed loading model" in out)
            }
        }
    }


    // sample random coverage and test the same model prediction at least.
    @Test
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
            val path = dir / "track.bed.gz"

            sampleCoverage(path, TO, BIN)
            print("Saved sampled track file: $path")

            val chromsizes = Genome["to1"].chromSizesPath.toString()
            SpanCLA.main(arrayOf("analyze",
                    "-cs", chromsizes,
                    "--workdir", dir.toString(),
                    "-t", path.toString(),
                    "--threads", THREADS.toString()))
            // Log file
            assertTrue((dir /  "logs" / "track_200.log").exists)
            // Coverage test
            assertTrue((dir / "cache" / "coverage").exists)
            assertEquals(1, (dir / "cache" / "coverage").glob("*").size)
            val coverageName = (dir / "cache" / "coverage").glob("*").first().fileName.toString()
            assertTrue("track_200_unique#.*\\.bw".toRegex().matches(coverageName))
            // Model test
            assertTrue((dir / "fit").exists)
            assertEquals(1, (dir / "fit").glob("*").size)
            val modelName = (dir / "fit").glob("*").first().fileName.toString()
            assertTrue("span_track_200_unique.tar".toRegex().matches(modelName))
            assertEquals(3, dir.glob("**/*.*").size)
        }

    }


    // sample random coverage and test the same model prediction at least.
    @Test
    fun testPeaksStats() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                SpanCLA.main(arrayOf("analyze",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString(),
                        "-o", path.toString()))
                val out = String(stream.toByteArray())
                assertFalse("""NO output path given, process model fitting only.
LABELS, FDR, GAP options are ignored.
""" in out)
                assertIn("Track source: $path", out)
                assertIn("Peaks Statistics:\n", out)
                assertIn("FRIP: ", out)
                // TODO[shpynov] fix me
                assertIn("Failed to estimate signal-to-noise ratio", out)
            }
        }
    }


    @Test
    fun analyzeSampledEnrichment() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) {
                val enriched = BitSet()
                if (it.name == "chr1") {
                    enriched.set(1000, 2000)
                }
                enriched
            }

            val zeroRegions = genomeMap(TO) {
                val zeroes = BitSet()
                if (it.name == "chr1") {
                    zeroes[3000] = 4000
                }
                zeroes
            }
            sampleCoverage(path, TO, BIN, enrichedRegions, zeroRegions)
            println("Saved sampled track file: $path")
            withTempDirectory("work") {
                val bedPath = it / "result.bed"
                SpanCLA.main(arrayOf("analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-only", "chr1",
                        "-w", it.toString(),
                        "-o", bedPath.toString(),
                        "-fdr", FDR.toString(),
                        "-t", path.toString()))
                SpanCLA.main(arrayOf("analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-only", "chr1",
                        "-w", it.toString(),
                        "-o", bedPath.toString(),
                        "-fdr", FDR.toString(),
                        "-t", path.toString()))
                // Check created bed file
                assertTrue(Location(1100 * BIN, 1900 * BIN, TO.get().first()) in LocationsMergingList.load(TO, bedPath))
                // Check correct log file name
                assertTrue((it / "logs" / "${bedPath.stem}.log").exists)
            }
        }
    }


    companion object {
        private val TO = GenomeQuery("to1")
        private val BIN = 200
        private val FDR = 1e-10
        private val GAP = 10
        private val THREADS = 1
        private val FRAGMENT = 150
        private val OUT = System.out
        private val ERR = System.err
        private val FORMAT = CSVFormat.TDF

        fun assertLinesEqual(expected: String, actual: String) = assertEquals(expected.lines(), actual.lines())

        fun assertIn(substring: String, fullString: String) {
            // Process Windows with different line separators correctly.
            for (s in substring.lines()) {
                assertTrue(s in fullString, "Expected <$s> to be in <$fullString>.")
            }
        }
    }
}
