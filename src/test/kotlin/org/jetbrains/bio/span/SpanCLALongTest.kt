package org.jetbrains.bio.span

import kotlinx.support.jdk7.use
import org.apache.commons.csv.CSVFormat
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
import org.jetbrains.bio.query.readsName
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.ByteArrayOutputStream
import java.io.FileReader
import java.io.PrintStream
import java.nio.file.Path
import java.security.Permission
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFails
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class SpanCLALongTest {

    @Before
    fun setUp() {
        (LogManager.getCurrentLoggers().toList() + listOf(LogManager.getRootLogger())).forEach {
            (it as Logger).level = Level.INFO
        }
        SpanCLA.ignoreConfigurePaths = true
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
                val bedPath = it / "peaks.bed"
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
                assertIn("""WARN Span] After fitting the model, emission's parameter p in LOW state
is higher than emission's parameter p in HIGH state
WARN Span] This is generally harmless, but could indicate low quality of data.
""", out)
            }
        }
    }


    @Test
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
            withTempFile("track", ".bed.gz", dir) { path ->

                sampleCoverage(path, TO, BIN, goodQuality = true)
                print("Saved sampled track file: $path")

                val chromsizes = Genome["to1"].chromSizesPath.toString()
                SpanCLA.main(arrayOf("analyze",
                        "-cs", chromsizes,
                        "--workdir", dir.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString()))

                // Check that log file was created correctly
                assertTrue((dir / "logs" / "${path.readsName()}_200.log").exists)

                /**
                 * Shpynov: since [Configuration] allows only single initialization of experimentPath,
                 * [SpanCLA] uses ignoreConfigurePaths variable to ignore it
                 */
                val dir = Configuration.experimentsPath

                // Coverage test
                assertTrue((dir / "cache" / "coverage").exists)
                assertTrue((dir / "cache" / "coverage").glob("{${path.readsName()}}_200_unique#*.bw").isNotEmpty())
                // Model test
                assertTrue((dir / "fit").exists)
                assertTrue((dir / "fit" / "${path.readsName()}_200_unique.span").exists)
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
                        "-o", peaksPath.toString()))
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

    @Test
    fun analyzeEmptyCoverage() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) { BitSet() }

            val zeroRegions = genomeMap(TO) {
                val zeroes = BitSet()
                zeroes.set(0, it.length / BIN)
                zeroes
            }

            sampleCoverage(path, TO, BIN, enrichedRegions, zeroRegions, goodQuality = true)
            println("Saved sampled track file: $path")

            val securityManager = object : SecurityManager() {
                override fun checkPermission(permission: Permission) {
                    if ("exitVM" in permission.name) {
                        throw SecurityException()
                    }
                }
            }

            withTempDirectory("work") {
                withSecurityManager(securityManager) {
                    try {
                        SpanCLA.main(arrayOf("analyze",
                                             "-cs", Genome["to1"].chromSizesPath.toString(),
                                             "-w", it.toString(),
                                             "-t", path.toString()))
                    } catch (ignore: Throwable) { }
                }

                // Check correct log file name
                val logPath = it / "logs" / "${reduceIds(listOf(path.readsName()))}_${BIN}.log"
                assertTrue(logPath.exists)
                assertIn("Model can't be trained on empty coverage, exiting.",
                         FileReader(logPath.toFile()).use { it.readText() })
            }
        }
    }


    companion object {
        private val TO = GenomeQuery("to1")
        private const val BIN = 200
        private const val FDR = 1e-10
        private const val GAP = 10
        private const val THREADS = 1
        private const val FRAGMENT = 150
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
                            while (true) {
                                val sample = model.sample(bins).sliceAsInt("d0")
                                if (sample.any { it != 0 }) return@run sample
                            }
                            return@run IntArray(0)
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

        inline fun withSecurityManager(securityManager: SecurityManager, block: () -> Any) {
            val oldManager = System.getSecurityManager()
            System.setSecurityManager(securityManager)
            block.invoke()
            System.setSecurityManager(oldManager)
        }

    }
}
