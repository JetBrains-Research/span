package org.jetbrains.bio.span

import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.span.coverage.SpanCoverageSampler.sampleCoverage
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BIN
import org.jetbrains.bio.span.fit.SpanModelType
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.text.DecimalFormatSymbols
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * RM - regression mixture
 */
class SpanRegrMixtureCLALongTest {

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
    fun testMultipleTracks() {
        withTempFile("trackA", ".bed.gz") { pathA ->
            sampleCoverage(pathA, SpanCLALongTest.TO, SPAN_DEFAULT_BIN, goodQuality = true)
            println("Saved sampled track file: $pathA")

            withTempFile("trackB", ".bed.gz") { pathB ->
                sampleCoverage(
                    pathA,
                    SpanCLALongTest.TO,
                    SPAN_DEFAULT_BIN,
                    goodQuality = false
                )
                println("Saved sampled track file: $pathA")

                withTempDirectory("work") { dir ->
                    val chromsizes = Genome["to1"].chromSizesPath.toString()

                    val (_, err) = Logs.captureLoggingOutput {
                        SpanCLA.main(
                            arrayOf(
                                "analyze",
                                "-cs", chromsizes,
                                "--workdir", dir.toString(),
                                "-t", listOf(pathA, pathB).joinToString(","),
                                "--threads", SpanCLALongTest.THREADS.toString(),
                                "--model-type", SpanModelType.POISSON_REGRESSION_MIXTURE.id
                            )
                        )
                    }
                    assertIn("ERROR", err)
                    assertIn("Poisson regression mixture currently accepts a single data track", err)
                }
            }
        }
    }

    @Test
    fun testOutput() {
        // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
        withTempFile("track", ".bed.gz") { path ->

            sampleCoverage(path, SpanCLALongTest.TO, SPAN_DEFAULT_BIN, goodQuality = true)
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
                            "--threads", SpanCLALongTest.THREADS.toString(),
                            "--model-type", SpanModelType.POISSON_REGRESSION_MIXTURE.id,
                            "--peaks", peaksPath.toString(),
                            "--fdr", "0.01"
                        )
                    )
                }
                assertFalse(
                    """NO output path given, process model fitting only.
    LABELS, FDR, GAP options are ignored.
    """ in out
                )
                val ds =
                    DecimalFormatSymbols(Locale.getDefault()).decimalSeparator // XXX: Not so important to make to types of tests for US and EU locales
                assertIn("100${ds}00% (", out)
                assertIn("MODEL TYPE: ${SpanModelType.POISSON_REGRESSION_MIXTURE}", out)
                assertIn("File: $peaksPath", out)
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
                assertIn("Signal to noise: ", out)
            }
        }
    }

    @Test
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    sampleCoverage(
                        path,
                        SpanCLALongTest.TO,
                        SPAN_DEFAULT_BIN,
                        goodQuality = true
                    )
                    sampleCoverage(
                        control,
                        SpanCLALongTest.TO,
                        SPAN_DEFAULT_BIN,
                        goodQuality = false
                    )

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--model-type", SpanModelType.POISSON_REGRESSION_MIXTURE.id,
                            "--threads", SpanCLALongTest.THREADS.toString(),
                            "--kc"
                        )
                    )

                    // Check that log file was created correctly
                    assertTrue(
                        (dir / "logs" / "${
                            reduceIds(
                                listOf(
                                    path.stemGz, control.stemGz,
                                    SPAN_DEFAULT_BIN.toString(), "unique"
                                )
                            )
                        }.log").exists,
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
                                "${reduceIds(listOf(path.stemGz, control.stemGz, SPAN_DEFAULT_BIN.toString()))}." +
                                        SpanModelType.POISSON_REGRESSION_MIXTURE.extension
                            ).size
                    )
                }
            }
        }
    }

    /**
     * Model extension is used to determine the model type.
     * @see [SpanModelType] for details
     * If the model extension contradicts the provided '--model-type' command line argument,
     * Span should exit with an error.
     */
    @Test
    fun testCustomModelPathWrongExtension() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                sampleCoverage(path, SpanCLALongTest.TO, SPAN_DEFAULT_BIN, goodQuality = true)

                val chromsizes = Genome["to1"].chromSizesPath.toString()

                val defaultModelPath = dir / "custom" / "path" / "model.span"
                val (_, wrongErr) = Logs.captureLoggingOutput {
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "--threads", SpanCLALongTest.THREADS.toString(),
                            "--model-type", SpanModelType.POISSON_REGRESSION_MIXTURE.id,
                            "--model", defaultModelPath.toString()
                        )
                    )
                }
                assertIn(
                    "model type (${SpanModelType.NB2Z_HMM}) " +
                            "differs from the command line argument (${SpanModelType.POISSON_REGRESSION_MIXTURE})",
                    wrongErr
                )
            }
        }
    }

}