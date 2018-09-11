package org.jetbrains.bio.span

import kotlinx.support.jdk7.use
import org.apache.log4j.Level
import org.apache.log4j.LogManager
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.Logs
import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.ByteArrayOutputStream
import java.io.FileReader
import java.io.PrintStream
import java.nio.file.Path
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
        SpanCLA.ignoreConfigurePaths = true
        Sampling.RANDOM_DATA_GENERATOR.randomGenerator.setSeed(1234L)
    }

    @After
    fun tearDown() {
        SpanCLA.ignoreConfigurePaths = false
        System.setOut(OUT)
        System.setErr(ERR)
        // we might have unfinished tracked tasks which will never be complete, let's drop them
        MultitaskProgress.clear()
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
compare                         Differential peak calling mode, experimental
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
compare                         Differential peak calling mode, experimental
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
compare                         Differential peak calling mode, experimental
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
compare                         Differential peak calling mode, experimental
""".trim(), String(stream.toByteArray()).trim())
    }


    @Test
    fun checkVersion() {
        val stream = ByteArrayOutputStream()
        System.setOut(PrintStream(stream))
        SpanCLA.main(arrayOf("--version"))
        assertEquals("@VERSION@.@BUILD@ built on @DATE@", String(stream.toByteArray()).trim())
    }


    @Test
    fun compareSameTestOrganismTracks() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = true)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val peaksPath = it / "peaks.bed"
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                // Update with changed System.out
                Logs.addConsoleAppender(Level.INFO)

                val chromsizes = Genome["to1"].chromSizesPath.toString()
                SpanCLA.main(arrayOf("compare",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t1", path.toString(),
                        "-t2", path.toString(),
                        "--peaks", peaksPath.toString(),
                        "--fdr", FDR.toString(),
                        "--gap", GAP.toString(),
                        "--threads", THREADS.toString()))

                assertTrue(peaksPath.size.isEmpty(), "Found differential peaks in identical signals.")

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
PEAKS: $peaksPath
""", out)
                assertIn("Saved result to $peaksPath", out)
                // Check model fit has a progress
                assertIn("] 0.00% (0/250), Elapsed time", out)
            }
        }
    }

    @Test
    fun compareSameTestOrganismTracksReplicates() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = true)
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
                        "--peaks", bedPath.toString(),
                        "--fdr", FDR.toString()))
                assertTrue(bedPath.size.isEmpty(),
                        "Found differential peaks in identical signals.")
            }
        }
    }

    @Test
    fun testModelFitFile() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = true)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                // Update with changed System.out
                Logs.addConsoleAppender(Level.INFO)

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
                assertIn(".span: done in ", out)
                assertIn("Model saved: ", out)
                assertFalse("Loading model" in out)
                assertFalse("Completed loading model" in out)
            }
        }
    }

    @Test
    fun testBadTrackQualityWarning() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = false)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                // Update with changed System.out
                Logs.addConsoleAppender(Level.INFO)

                val chromsizes = Genome["to1"].chromSizesPath.toString()
                SpanCLA.main(arrayOf("analyze",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString()))
                val out = String(stream.toByteArray())
                assertIn("WARN SPAN] After fitting the model, emission's parameter p in LOW state", out)
                assertIn("WARN SPAN] This is generally harmless, but could indicate low quality of data.", out)
            }
        }
    }


    @Test
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(path, TO, BIN, goodQuality = true)
                    sampleCoverage(control, TO, BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    SpanCLA.main(arrayOf("analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString()))

                    // Check that log file was created correctly
                    assertTrue((dir / "logs" / "${path.stemGz}_${control.stemGz}_200.log").exists)

                    assertTrue((Configuration.experimentsPath / "cache").exists)

                    // Genome Coverage test
                    assertEquals(1, (Configuration.experimentsPath / "cache").glob("coverage_${path.stemGz}_unique#*.npz").size)
                    assertEquals(1, (Configuration.experimentsPath / "cache").glob("coverage_${control.stemGz}_unique#*.npz").size)
                    // Model test
                    assertTrue((Configuration.experimentsPath / "fit").exists)
                    assertEquals(1, (Configuration.experimentsPath / "fit").glob("${path.stemGz}_${control.stemGz}_200#*.span").size)
                }
            }
        }
    }


    @Test
    fun testPeaksStats() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val stream = ByteArrayOutputStream()
                System.setOut(PrintStream(stream))
                // Update with changed System.out
                Logs.addConsoleAppender(Level.INFO)

                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val peaksPath = path.parent / "${path.stem}.peak"
                SpanCLA.main(arrayOf("analyze", "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString(),
                        "--peaks", peaksPath.toString()))
                val out = String(stream.toByteArray())
                assertFalse("""NO output path given, process model fitting only.
    LABELS, FDR, GAP options are ignored.
    """ in out)
                assertIn("Track source: $peaksPath", out)
                assertIn("Peaks Statistics:\n", out)
                assertIn("FRIP: ", out)
                assertIn("Signal to noise: ", out)
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
            sampleCoverage(path, TO, BIN, enrichedRegions, zeroRegions, goodQuality = true)
            println("Saved sampled track file: $path")
            withTempDirectory("work") {
                val bedPath = it / "result.bed"
                SpanCLA.main(arrayOf("analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", it.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", FDR.toString(),
                        "-t", path.toString()))
                SpanCLA.main(arrayOf("analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", it.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", FDR.toString(),
                        "-t", path.toString()))
                // Check created bed file
                assertTrue(Location(1100 * BIN, 1900 * BIN, TO.get().first()) in LocationsMergingList.load(TO, bedPath))
                // Check correct log file name
                assertTrue((it / "logs" / "${bedPath.stem}.log").exists)
            }
        }
    }

    @Test
    fun analyzeEmptyCoverage() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) { BitSet() }

            val zeroRegions = genomeMap(TO) {
                val zeroes = BitSet()
                zeroes.set(0, it.length / BIN)
                zeroes
            }

            val outStream = ByteArrayOutputStream()
            val errStream = ByteArrayOutputStream()
            System.setOut(PrintStream(outStream))
            System.setErr(PrintStream(errStream))
            // Update with changed System.out
            Logs.addConsoleAppender(Level.INFO)

            sampleCoverage(path, TO, BIN, enrichedRegions, zeroRegions, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                /* Turn suppressExit on, otherwise Span would call System.exit */
                withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                    SpanCLA.main(arrayOf("analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", it.toString(),
                            "-t", path.toString()))
                }

                val out = String(outStream.toByteArray())
                val err = String(errStream.toByteArray())

                // Check correct log file name
                val logPath = it / "logs" / "${reduceIds(listOf(path.stemGz))}_${BIN}.log"
                assertTrue(logPath.exists)
                val log = FileReader(logPath.toFile()).use { it.readText() }
                val errorMessage = "Model can't be trained on empty coverage, exiting."
                assertIn(errorMessage, log)
                assertIn(errorMessage, out)
                assertIn(errorMessage, err)
            }
        }
    }

    companion object {
        private val TO = GenomeQuery("to1")
        private const val BIN = 200
        private const val FDR = 1E-10
        private const val GAP = 10
        private const val THREADS = 1
        private const val FRAGMENT = 150
        private val OUT = System.out
        private val ERR = System.err

        fun assertLinesEqual(expected: String, actual: String) = assertEquals(expected.lines(), actual.lines())

        fun assertIn(substring: String, fullString: String) {
            // Process Windows with different line separators correctly.
            for (s in substring.lines()) {
                assertTrue(s in fullString, "Expected <$s> to be in <$fullString>.")
            }
        }


        fun sampleCoverage(path: Path, genomeQuery: GenomeQuery, bin: Int, goodQuality: Boolean) =
                sampleCoverage(path,
                        genomeQuery, bin,
                        genomeMap(genomeQuery) { BitSet() },
                        genomeMap(genomeQuery) { BitSet() },
                        goodQuality)

        fun sampleCoverage(path: Path,
                           genomeQuery: GenomeQuery, bin: Int,
                           fulls: GenomeMap<BitSet>, zeroes: GenomeMap<BitSet>,
                           goodQuality: Boolean) {
            withResource(SpanCLALongTest::class.java,
                    if (goodQuality)
                        "GSM646345_H1_H3K4me3_rep1_hg19_model.json"
                    else
                        "yd6_k27ac_failed_model.json") { modelPath ->
                val model = ClassificationModel.load<MLFreeNBHMM>(modelPath)
                BedFormat().print(path).use {
                    genomeQuery.get().forEach { chr ->
                        val bins = chr.length / bin
                        val coverage = run {
                            (0..9).forEach {
                                val sample = model.sample(bins).sliceAsInt("d0")
                                if (sample.any { it != 0 }) return@run sample
                            }
                            throw IllegalStateException("Couldn't sample coverage.")
                        }
                        for (b in 0 until bins) {
                            val enriched = fulls[chr][b]
                            val zero = zeroes[chr][b]
                            check(!(enriched && zero)) { "Both enriched and zero for $b" }
                            val c = if (enriched) bin else if (zero) 0 else coverage[b]
                            for (i in 0 until c) {
                                val start = b * bin + i
                                it.print(BedEntry(chr.name, start, start + 1))
                            }
                        }
                    }
                }
            }
        }

        inline fun withSystemProperty(property: String, value: String, block: () -> Any) {
            val oldValue = System.setProperty(property, value)
            try {
                block.invoke()
            } finally {
                if (oldValue != null) System.setProperty(property, oldValue) else System.clearProperty(property)
            }
        }
    }
}
