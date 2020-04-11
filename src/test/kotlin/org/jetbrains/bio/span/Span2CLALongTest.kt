package org.jetbrains.bio.span

import org.jetbrains.bio.Tests
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.query.stemGz
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

class Span2CLALongTest {

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
    fun testMultipleTracks() {
        withTempFile("trackA", ".bed.gz") { pathA ->
            SpanCLALongTest.sampleCoverage(pathA, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = true)
            println("Saved sampled track file: $pathA")

            withTempFile("trackB", ".bed.gz") { pathB ->
                SpanCLALongTest.sampleCoverage(pathA, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = false)
                println("Saved sampled track file: $pathA")

                withTempDirectory("work") { dir ->
                    val chromsizes = Genome["to1"].chromSizesPath.toString()

                    val (_, err) = Logs.captureLoggingOutput {
                        SpanCLALongTest.withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                            SpanCLA.main(arrayOf(
                                "analyze-experimental",
                                "-cs", chromsizes,
                                "--workdir", dir.toString(),
                                "-t", listOf(pathA, pathB).joinToString(","),
                                "--threads", SpanCLALongTest.THREADS.toString(),
                                "--type", "prm"
                            ))
                        }
                    }
                    Tests.assertIn("ERROR", err)
                    Tests.assertIn("Poisson regression mixture currently accepts a single data track", err)
                }
            }
        }
    }

    @Test
    fun testOutput() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            SpanCLALongTest.sampleCoverage(path, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = true)
            println("Saved sampled track file: $path")

            withTempDirectory("work") {
                val chromsizes = Genome["to1"].chromSizesPath.toString()
                val peaksPath = path.parent / "${path.stem}.peak"
                val (out, _) = Logs.captureLoggingOutput {
                    SpanCLA.main(arrayOf("analyze-experimental", "-cs", chromsizes,
                        "--workdir", it.toString(),
                        "-t", path.toString(),
                        "--threads", SpanCLALongTest.THREADS.toString(),
                        "--type", "prm",
                        "--peaks", peaksPath.toString(),
                        "--fdr", "0.01"))
                }
                assertFalse("""NO output path given, process model fitting only.
    LABELS, FDR, GAP options are ignored.
    """ in out)
                Tests.assertIn("100.00% (", out)
                Tests.assertIn("MODEL TYPE: ${SpanModel.POISSON_REGRESSION_MIXTURE}", out)
                Tests.assertIn("Source: $peaksPath", out)
                Tests.assertIn("FRIP: ", out)

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
                Tests.assertIn("Signal to noise: ", out)
            }
        }
    }

    @Test
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    SpanCLALongTest.sampleCoverage(path, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = true)
                    SpanCLALongTest.sampleCoverage(control, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    SpanCLA.main(arrayOf(
                        "analyze-experimental",
                        "-cs", chromsizes,
                        "--workdir", dir.toString(),
                        "-t", path.toString(),
                        "-c", control.toString(),
                        "--type", "prm",
                        "--threads", SpanCLALongTest.THREADS.toString()
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
                                .glob("${reduceIds(listOf(path.stemGz, control.stemGz, "200"))}.span2").size
                    )
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
                SpanCLALongTest.sampleCoverage(path, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = true)

                val chromsizes = Genome["to1"].chromSizesPath.toString()

                val wrongModelPath = dir / "custom" / "path" / "model.span"
                val (_, wrongErr) = Logs.captureLoggingOutput {
                    SpanCLALongTest.withSystemProperty(JOPTSIMPLE_SUPPRESS_EXIT, "true") {
                        SpanCLA.main(arrayOf(
                            "analyze-experimental",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", SpanCLALongTest.THREADS.toString(),
                            "--type", "prm",
                            "--model", wrongModelPath.toString()
                        ))
                    }
                }
                Tests.assertIn(
                    "Stored model type (${SpanModel.NB_HMM}) " +
                            "differs from the command line argument (${SpanModel.POISSON_REGRESSION_MIXTURE})",
                    wrongErr
                )
            }
        }
    }

}