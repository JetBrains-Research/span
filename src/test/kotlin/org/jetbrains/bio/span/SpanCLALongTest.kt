package org.jetbrains.bio.span

import kotlinx.support.jdk7.use
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.big.BedEntry
import org.jetbrains.bio.experiments.fit.SpanModel
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
import java.io.FileReader
import java.nio.file.Path
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertNotEquals
import kotlin.test.assertTrue

class SpanCLALongTest {

    @Before
    fun setUp() {
        SpanCLA.ignoreConfigurePaths = true
        Sampling.RANDOM_DATA_GENERATOR.randomGenerator.setSeed(1234L)
    }

    @After
    fun tearDown() {
        SpanCLA.ignoreConfigurePaths = false
        // we might have unfinished tracked tasks which will never be complete, let's drop them
        MultitaskProgress.clear()
    }

    @Test
    fun emptyArgs() {
        val (_, err) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf())
        }
        assertLinesEqual("""
ERROR: No command given; analyze or compare expected.

Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode, experimental
""".trim(), err.trim())
    }

    @Test
    fun illegalArgs() {
        val (_, err) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("foobar"))
        }
        assertLinesEqual("""
ERROR: Unknown command: foobar; analyze or compare expected.

Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode, experimental
""".trim(), err.trim())
    }

    @Test
    fun quietError() {
        val (_, err) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("foobar", "quiet"))
        }
        assertLinesEqual("""
ERROR: Unknown command: foobar; analyze or compare expected.

Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode, experimental
        """.trim(), err.trim())
    }

    @Test
    fun checkHelp() {
        val (out, _) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("--help"))
        }
        assertLinesEqual("""
Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode, experimental
""".trim(), out.trim())
    }


    @Test
    fun checkVersion() {
        val (out, _) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("--version"))
        }
        assertEquals("@VERSION@.@BUILD@ built on @DATE@", out.trim())
    }


    @Test
    fun compareSameTestOrganismTracks() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = true)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val peaksPath = it / "peaks.bed"
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(arrayOf(
                        "compare",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t1", path.toString(),
                        "-t2", path.toString(),
                        "--peaks", peaksPath.toString(),
                        "--fdr", FDR.toString(),
                        "--gap", GAP.toString(),
                        "--threads", THREADS.toString()
                    ))
                }

                assertTrue(
                    peaksPath.size.isEmpty(),
                    "Found differential peaks in identical signals."
                )

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
FRAGMENT: auto
BIN: $BIN
FDR: $FDR
GAP: $GAP
PEAKS: $peaksPath
""", out)
                assertIn("Saved result to $peaksPath", out)
                // Check model fit has a progress
                assertIn("] 0.00% (0/100), Elapsed time", out)
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
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(arrayOf(
                        "analyze",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString()
                    ))
                }
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
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(arrayOf(
                        "analyze",
                        "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString()
                    ))
                }
                assertIn(
                    "] WARN Span After fitting the model, emission's parameter p in LOW state", out
                )
                assertIn(
                    "] WARN Span This is generally harmless, but could indicate low quality of data.", out
                )
            }
        }
    }

    @Test
    fun testQuietMode() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = false)
            print("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                /* We can't test the stdout, because the Span quiet mode redefines the JVM [System.out].
                 * But we can restore the System.out to the original value using [captureLoggingOutput].
                 */
                Logs.captureLoggingOutput {
                    val oldSystemOut = System.out
                    SpanCLA.main(arrayOf(
                        "analyze",
                        "-cs", chromsizes,
                        "--workdir", dir.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString(),
                        "-q"
                    ))
                    /* the best we can do is to check whether any redirection took place */
                    assertNotEquals(
                        oldSystemOut, System.out,
                        "Span quiet mode didn't redirect System.out"
                    )
                }
                /* we also check that logging was performed normally */
                val logPath = dir / "logs" / "${reduceIds(listOf(path.stemGz, BIN.toString(), "unique"))}.log"
                assertTrue(logPath.exists, "Log file not found")
                assertTrue(logPath.size.isNotEmpty(), "Log file is empty")
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
                    SpanCLA.main(arrayOf(
                        "analyze",
                        "-cs", chromsizes,
                        "--workdir", dir.toString(),
                        "-t", path.toString(),
                        "-c", control.toString(),
                        "--threads", THREADS.toString()
                    ))

                    // Check that log file was created correctly
                    assertTrue(
                        (dir / "logs" / "${reduceIds(listOf(path.stemGz, control.stemGz, "200", "unique"))}.log")
                                .exists,
                        "Log file not found"
                    )

                    assertTrue((Configuration.experimentsPath / "cache").exists)

                    // Genome Coverage test
                    assertEquals(
                        1,
                        (Configuration.experimentsPath / "cache")
                                .glob("coverage_${path.stemGz}_unique#*.npz").size
                    )
                    assertEquals(
                        1,
                        (Configuration.experimentsPath / "cache")
                                .glob("coverage_${control.stemGz}_unique#*.npz").size)
                    // Model test
                    assertTrue((Configuration.experimentsPath / "fit").exists)
                    assertEquals(
                        1,
                        (Configuration.experimentsPath / "fit")
                                .glob("${reduceIds(listOf(path.stemGz, control.stemGz, "200"))}.span").size
                    )
                }
            }
        }
    }

    @Test
    fun testCustomModelPath() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(path, TO, BIN, goodQuality = true)
                    sampleCoverage(control, TO, BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    val modelPath = dir / "custom" / "path" / "model.span"
                    SpanCLA.main(arrayOf(
                        "analyze",
                        "-cs", chromsizes,
                        "--workdir", dir.toString(),
                        "-t", path.toString(),
                        "-c", control.toString(),
                        "--threads", THREADS.toString(),
                        "--model", modelPath.toString()
                    ))
                    assertTrue(modelPath.exists, "Model was not created at $modelPath")
                    assertTrue(modelPath.size.isNotEmpty(), "Model file $modelPath is empty")
                    val (reloadOut, reloadErr) = Logs.captureLoggingOutput {
                        SpanCLA.main(arrayOf(
                            "analyze",
                            "--workdir", dir.toString(),
                            "--threads", THREADS.toString(),
                            "--model", modelPath.toString()
                        ))
                    }
                    assertIn("Completed loading model: $modelPath", reloadOut)
                    assertEquals("", reloadErr)

                    val (_, invalidErr) = Logs.captureLoggingOutput {
                        withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                            SpanCLA.main(arrayOf("analyze",
                                "--workdir", dir.toString(),
                                "--threads", THREADS.toString(),
                                "--model", modelPath.toString(),
                                "--bin", "137"
                            ))
                        }
                    }
                    assertIn("Stored bin size (200) differs from the command line argument (137)", invalidErr)
                }
            }
        }
    }

    /**
     * Model extension is used to determine the model type.
     * .span = negative binomial HMM (classic Span)
     * .span2 = Poisson regression mixture (experimental Span)
     * any other = error, unrecognized type.
     * If the model extension contradicts the provided '--type' command line argument, Span should exit with an error.
     */
    @Test
    fun testCustomModelPathWrongExtension() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, TO, BIN, goodQuality = true)

                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val invalidModelPath = dir / "custom" / "path" / "model.foo"
                val (_, invalidErr) = Logs.captureLoggingOutput {
                    withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                        SpanCLA.main(arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--model", invalidModelPath.toString()
                        ))
                    }
                }
                assertIn(
                    "Unrecognized model extension '.foo', should be either '.span' or '.span2'",
                    invalidErr
                )

                val wrongModelPath = dir / "custom" / "path" / "model.span2"
                val (_, wrongErr) = Logs.captureLoggingOutput {
                    withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                        SpanCLA.main(arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--type", "nbhmm",
                            "--model", wrongModelPath.toString()
                        ))
                    }
                }
                assertIn(
                    "Stored model type (${SpanModel.POISSON_REGRESSION_MIXTURE}) " +
                        "differs from the command line argument (${SpanModel.NB_HMM})",
                    wrongErr
                )
            }
        }
    }


    @Test
    fun testOutput() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val peaksPath = path.parent / "${path.stem}.peak"
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(arrayOf("analyze", "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", THREADS.toString(),
                        "--peaks", peaksPath.toString()))
                }
                assertFalse("""NO output path given, process model fitting only.
    LABELS, FDR, GAP options are ignored.
    """ in out)
                assertIn("100.00% (", out)
                assertIn("Track source: $peaksPath", out)
                assertIn("FRIP: ", out)

                /* Check that coverage is being generated */
                val format = BedFormat.from("bed6+3")
                assertTrue(
                    format.parse(peaksPath) { parser ->
                        parser.all { entry ->
                            val coverage = entry.unpack(6).extraFields?.get(0)
                            return@all coverage != null && coverage != "0.0"
                        }
                    },
                    "Peak value is reported as 0.0, although the coverage cache is present"
                )
                assertIn("Multistart done", out)
                assertIn("Signal mean: ", out)
                assertIn("Noise mean: ", out)
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
            withTempDirectory("work") { dir ->
                val bedPath = dir / "result.bed"
                SpanCLA.main(arrayOf(
                    "analyze",
                    "-cs", Genome["to1"].chromSizesPath.toString(),
                    "-w", dir.toString(),
                    "--peaks", bedPath.toString(),
                    "-fdr", FDR.toString(),
                    "-t", path.toString()
                ))
                SpanCLA.main(arrayOf(
                    "analyze",
                    "-cs", Genome["to1"].chromSizesPath.toString(),
                    "-w", dir.toString(),
                    "--peaks", bedPath.toString(),
                    "-fdr", FDR.toString(),
                    "-t", path.toString()
                ))
                // Check created bed file
                assertTrue(
                    Location(1100 * BIN, 1900 * BIN, TO.get().first())
                            in LocationsMergingList.load(TO, bedPath),
                    "Expected location not found in called peaks"
                )
                // Check correct log file name
                val logPath = dir / "logs" / "${bedPath.stem}.log"
                assertTrue(logPath.exists, "Log file not found")
                val log = FileReader(logPath.toFile()).use { it.readText() }
                assertIn("Signal mean:", log)
                assertIn("Noise mean:", log)
                assertIn("Signal to noise:", log)
                /** In sampling producer model - signal to noise ratio ~ 30
                 * "mean": 0.36580807929296383, // noise
                 * "mean": 9.113896687733767,   // signal
                 */
                assertTrue(log.substringAfter("Signal to noise:").substringBefore("\n").trim().toDouble() > 10)
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

            withTempDirectory("work") { dir ->
                /* Turn suppressExit on, otherwise Span would call System.exit */
                val (out, err) = Logs.captureLoggingOutput {
                    withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                        SpanCLA.main(arrayOf("analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "-t", path.toString()))
                    }
                }

                // Check correct log file name
                val logPath = dir / "logs" / "${reduceIds(listOf(path.stemGz, BIN.toString(), "unique"))}.log"
                assertTrue(logPath.exists, "Log file not found")
                val log = FileReader(logPath.toFile()).use { it.readText() }
                val errorMessage = "Model can't be trained on empty coverage, exiting."
                assertIn(errorMessage, log)
                assertIn(errorMessage, out)
                assertIn(errorMessage, err)
            }
        }
    }

    companion object {
        internal val TO = GenomeQuery(Genome["to1"])
        internal const val BIN = 200
        private const val FDR = 1E-10
        private const val GAP = 10
        internal const val THREADS = 1
        private const val FRAGMENT = 150

        fun assertLinesEqual(expected: String, actual: String) =
                assertEquals(expected.lines(), actual.lines())

        fun sampleCoverage(path: Path, genomeQuery: GenomeQuery, bin: Int, goodQuality: Boolean) =
                sampleCoverage(
                    path,
                    genomeQuery, bin,
                    genomeMap(genomeQuery) { BitSet() },
                    genomeMap(genomeQuery) { BitSet() },
                    goodQuality
                )

        fun sampleCoverage(
                path: Path,
                genomeQuery: GenomeQuery, bin: Int,
                fulls: GenomeMap<BitSet>, zeroes: GenomeMap<BitSet>,
                goodQuality: Boolean
        ) {
            withResource(SpanCLALongTest::class.java,
                if (goodQuality)
                    "GSM646345_H1_H3K4me3_rep1_hg19_model.json"
                else
                    "yd6_k27ac_failed_model.json") { modelPath ->
                val model = ClassificationModel.load<MLFreeNBHMM>(modelPath)
                model.logPriorProbabilities[0] = Double.NEGATIVE_INFINITY
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
