package org.jetbrains.bio.span

import kotlinx.support.jdk7.use
import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.Tests.assertMatches
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.span.coverage.SpanCoverageSampler.sampleCoverage
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BIN
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FDR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FIT_THRESHOLD
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.FileReader
import java.text.DecimalFormatSymbols
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
        System.setProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true")
    }

    @After
    fun tearDown() {
        SpanCLA.ignoreConfigurePaths = false
        // we might have unfinished tracked tasks which will never be complete, let's drop them
        MultitaskProgress.clear()
        System.setProperty(JOPTSIMPLE_SUPPRESS_EXIT, "false")
    }

    @Test
    fun emptyArgs() {
        val (_, err) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf())
        }
        assertTrue("ERROR: No command given." in err)
    }

    @Test
    fun illegalArgs() {
        val (_, err) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("foobar"))
        }
        assertTrue("ERROR: Unknown command: foobar." in err)
    }

    @Test
    fun quietError() {
        val (_, err) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("foobar", "quiet"))
        }
        assertTrue("ERROR: Unknown command: foobar." in err)
    }


    @Test
    fun checkVersion() {
        val (out, _) = Logs.captureLoggingOutput {
            SpanCLA.main(arrayOf("--version"))
        }
        // the test is sometimes launched in the assembled JAR, where the @@ tokens have already been substituted
        assertMatches(
            out.trim(),
            Regex("^[0-9]+\\.[0-9]+(\\.dev)?\\.[0-9]+ built on [A-Z][a-z]* [0-9]{2}, [0-9]{4}|@VERSION@\\.@BUILD@ built on @DATE@$")
        )
    }


    @Test
    fun compareSameTestOrganismTracks() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val peaksPath = it / "peaks.bed"
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "compare",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t1", path.toString(),
                            "-t2", path.toString(),
                            "--peaks", peaksPath.toString(),
                            "--fdr", SPAN_DEFAULT_FDR.toString(),
                            "--gap", SPAN_DEFAULT_GAP.toString(),
                            "--threads", THREADS.toString()
                        )
                    )
                }

                assertTrue(
                    peaksPath.size.isEmpty(),
                    "Found differential peaks in identical signals."
                )

                assertIn(
                    """SPAN
COMMAND:
LOG:
WORKING DIR: $it
THREADS: $THREADS
TREATMENT1: $path
CONTROL1: none
TREATMENT2: $path
CONTROL2: none
CHROM.SIZES: $chromsizes
FRAGMENT: auto
BIN: $SPAN_DEFAULT_BIN
FDR: $SPAN_DEFAULT_FDR
GAP: $SPAN_DEFAULT_GAP
PEAKS: $peaksPath
""", out
                )
                assertIn("Saved result to $peaksPath", out)

                // Check model fit has progress:
                // XXX: Not so important to make to types of tests for US and EU locales
                val ds = DecimalFormatSymbols(Locale.getDefault()).decimalSeparator
                assertIn("0${ds}00% (0/${SPAN_FIT_MAX_ITERATIONS}), Elapsed time", out)
                assertIn("100${ds}00% (", out)
            }
        }
    }

    @Test
    fun compareSameTestOrganismTracksReplicates() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val bedPath = it / "peaks.bed"
                SpanCLA.main(
                    arrayOf(
                        "compare",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", it.toString(),
                        "-b", SPAN_DEFAULT_BIN.toString(),
                        "-g", SPAN_DEFAULT_GAP.toString(),
                        "-fragment", FRAGMENT.toString(),
                        "-t1", "$path,$path",
                        "-t2", "$path,$path,$path",
                        "--peaks", bedPath.toString(),
                        "--fdr", SPAN_DEFAULT_FDR.toString()
                    )
                )
                assertTrue(
                    bedPath.size.isEmpty(),
                    "Found differential peaks in identical signals."
                )
            }
        }
    }

    @Test
    fun testModelFitFile() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString()
                        )
                    )
                }
                assertIn(
                    """NO peaks path given, process model fitting only.
LABELS, FDR, GAP options are ignored.
""", out
                )
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

            sampleCoverage(path, TO, 200, goodQuality = false)
            print("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--bin", "200",
                            "--it", "20",
                            "--debug"
                        )
                    )
                }
                assertTrue("After fitting the model, emission's parameter p in LOW state" in out)
                assertTrue("Low quality of data detected after fitting the model." in out)
            }
        }
    }

    @Test
    fun testQuietMode() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = false)
            print("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                /* We can't test the stdout, because the Span quiet mode redefines the JVM [System.out].
                 * But we can restore the System.out to the original value using [captureLoggingOutput].
                 */
                Logs.captureLoggingOutput {
                    val oldSystemOut = System.out
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "-q"
                        )
                    )
                    /* the best we can do is to check whether any redirection took place */
                    assertNotEquals(
                        oldSystemOut, System.out,
                        "Span quiet mode didn't redirect System.out"
                    )
                }
                /* we also check that logging was performed normally */
                val logPath = dir / "logs" / "${
                    reduceIds(
                        listOf(
                            path.stemGz,
                            SPAN_DEFAULT_BIN.toString(),
                            "unique"
                        )
                    )
                }.log"
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
                    sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, SPAN_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString()
                        )
                    )

                    // Check that log file was created correctly
                    assertTrue(
                        (dir / "logs" / "${
                            reduceIds(
                                listOf(
                                    path.stemGz,
                                    control.stemGz,
                                    SPAN_DEFAULT_BIN.toString(),
                                    "unique"
                                )
                            )
                        }.log")
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
                            .glob("coverage_${control.stemGz}_unique#*.npz").size
                    )
                    // Model test
                    assertTrue((Configuration.experimentsPath / "fit").exists)
                    assertEquals(
                        1,
                        (Configuration.experimentsPath / "fit")
                            .glob(
                                "${
                                    reduceIds(
                                        listOf(
                                            path.stemGz,
                                            control.stemGz,
                                            SPAN_DEFAULT_BIN.toString()
                                        )
                                    )
                                }.span"
                            ).size
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
                    sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, SPAN_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    val modelPath = dir / "custom" / "path" / "model.span"
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--model", modelPath.toString()
                        )
                    )
                    assertTrue(modelPath.exists, "Model was not created at $modelPath")
                    assertTrue(modelPath.size.isNotEmpty(), "Model file $modelPath is empty")
                    val (reloadOut, reloadErr) = Logs.captureLoggingOutput {
                        SpanCLA.main(
                            arrayOf(
                                "analyze",
                                "--workdir", dir.toString(),
                                "--threads", THREADS.toString(),
                                "--model", modelPath.toString()
                            )
                        )
                    }
                    assertIn("Completed loading model: ${modelPath.stem}", reloadOut)
                    assertEquals("", reloadErr)

                    val (_, invalidErr) = Logs.captureLoggingOutput {
                        SpanCLA.main(
                            arrayOf(
                                "analyze",
                                "--workdir", dir.toString(),
                                "--threads", THREADS.toString(),
                                "--model", modelPath.toString(),
                                "--bin", "137"
                            )
                        )
                    }
                    assertIn(
                        "bin size (${SPAN_DEFAULT_BIN}) differs from the command line argument (137)",
                        invalidErr
                    )
                }
            }
        }
    }

    /**
     * Classical Span only recognizes '.span' model file extension.
     * Other extensions are reserved for future use.
     */
    @Test
    fun testCustomModelPathWrongExtension() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)

                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val invalidModelPath = dir / "custom" / "path" / "model.foo"
                val (_, invalidErr) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--model", invalidModelPath.toString()
                        )
                    )
                }
            }
        }
    }

    @Test
    fun testTypeOnlyValidInExperimentalMode() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)

                val chromsizes = Genome["to1"].chromSizesPath.toString()

                val (_, wrongErr) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--type", "nbhmm"
                        )
                    )
                }
                assertIn(
                    "ERROR: type is not a recognized option",
                    wrongErr
                )
            }
        }
    }


    @Test
    fun testAnalyze() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, TO, SPAN_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val peaksPath = path.parent / "${path.stem}.peak"
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze", "-cs", chromsizes,
                            "--workdir", it.toString(),
                            "-t", path.toString(),
                            "--threads", THREADS.toString(),
                            "--peaks", peaksPath.toString()
                        )
                    )
                }

                val ds =
                    DecimalFormatSymbols(Locale.getDefault()).decimalSeparator // XXX: Not so important to make to types of tests for US and EU locales
                assertIn(
                    """SPAN
COMMAND:
LOG:
WORKING DIR: $it
THREADS: $THREADS
TREATMENT: $path
CONTROL: none
CHROM.SIZES: $chromsizes
FRAGMENT: auto
MAX ITERATIONS: $SPAN_FIT_MAX_ITERATIONS
CONVERGENCE THRESHOLD: $SPAN_FIT_THRESHOLD
CLIP: false
EXTENDED MODEL INFO: false
Library: ${path.fileName}, Depth:
100${ds}00% (
File: $peaksPath
FRIP: 
Reads: single-ended, Fragment size: 2 bp (cross-correlation estimate)
""", out
                )
                assertFalse(
                    """NO peaks path given, process model fitting only.
    LABELS, FDR, GAP options are ignored.
    """ in out
                )

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
            sampleCoverage(
                path,
                TO,
                SPAN_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            withTempDirectory("work") { dir ->
                val bedPath = dir / "result.bed"
                SpanCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", SPAN_DEFAULT_FDR.toString(),
                        "-t", path.toString()
                    )
                )
                // Check created bed file
                assertTrue(
                    Location(
                        1100 * SPAN_DEFAULT_BIN,
                        1900 * SPAN_DEFAULT_BIN,
                        TO.get().first()
                    )
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
    fun analyzeSampledEnrichmentPeaksVsIslands() {
        withTempFile("track", ".bed.gz") { coveragePath ->
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
            sampleCoverage(
                coveragePath,
                TO,
                SPAN_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $coveragePath")

            withTempDirectory("work") { dir ->
                val peaksPath = dir / "peaks.bed"
                SpanCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "--peaks", peaksPath.toString(),
                        "-fdr", SPAN_DEFAULT_FDR.toString(),
                        "-t", coveragePath.toString()
                    )
                )
                // Check created bed file
                val peaksLocations = LocationsMergingList.load(TO, peaksPath)
                assertTrue(
                    Location(
                        1100 * SPAN_DEFAULT_BIN,
                        1900 * SPAN_DEFAULT_BIN,
                        TO.get().first()
                    ) in peaksLocations,
                    "Expected location not found in called peaks"
                )
            }
        }
    }


    @Test
    fun analyzeSampledEnrichmentReusingModel() {
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
            sampleCoverage(
                path,
                TO,
                SPAN_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")
            withTempDirectory("work") { dir ->
                val bedPath = dir / "result.bed"
                val modelPath = dir / "model.span"
                SpanCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "-m", modelPath.toString(),
                        "-t", path.toString()
                    )
                )
                SpanCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-m", modelPath.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", SPAN_DEFAULT_FDR.toString()
                    )
                )
                // Check created bed file
                assertTrue(
                    Location(
                        1100 * SPAN_DEFAULT_BIN,
                        1900 * SPAN_DEFAULT_BIN,
                        TO.get().first()
                    ) in LocationsMergingList.load(TO, bedPath),
                    "Expected location not found in called peaks"
                )
            }
        }
    }


    @Test
    fun analyzeEmptyCoverage() {
        withTempFile("track", ".bed.gz") { path ->
            val enrichedRegions = genomeMap(TO) { BitSet() }

            val zeroRegions = genomeMap(TO) {
                val zeroes = BitSet()
                zeroes.set(0, it.length / SPAN_DEFAULT_BIN)
                zeroes
            }



            sampleCoverage(
                path,
                TO,
                SPAN_DEFAULT_BIN,
                enrichedRegions,
                zeroRegions,
                goodQuality = true
            )
            println("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                /* Turn suppressExit on, otherwise Span would call System.exit */
                val (out, err) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", Genome["to1"].chromSizesPath.toString(),
                            "-w", dir.toString(),
                            "-t", path.toString()
                        )
                    )
                }

                // Check correct log file name
                val logPath = dir / "logs" / "${
                    reduceIds(
                        listOf(
                            path.stemGz,
                            SPAN_DEFAULT_BIN.toString(),
                            "unique"
                        )
                    )
                }.log"
                assertTrue(logPath.exists, "Log file not found")
                val log = FileReader(logPath.toFile()).use { it.readText() }
                val errorMessage = "Model can't be trained on empty coverage, exiting."
                assertIn(errorMessage, log)
                assertIn(errorMessage, out)
                assertIn(errorMessage, err)
            }
        }
    }

    @Test
    fun analyzePartiallyEmptyCoverage() {
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(
                path,
                GenomeQuery(Genome["to1"], "chr1", "chr2"),
                SPAN_DEFAULT_BIN,
                goodQuality = true
            )
            println("Saved sampled track file: $path")

            withTempDirectory("work") { dir ->
                SpanCLA.main(
                    arrayOf(
                        "analyze",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "-t", path.toString()
                    )
                )
            }
        }
    }

    companion object {
        internal val TO = GenomeQuery(Genome["to1"])
        internal const val THREADS = 1
        private const val FRAGMENT = 200
    }
}
