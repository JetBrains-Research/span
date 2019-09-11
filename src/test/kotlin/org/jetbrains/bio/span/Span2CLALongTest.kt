package org.jetbrains.bio.span

import org.jetbrains.bio.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import kotlin.test.assertEquals
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
    fun testFilesCreatedByAnalyze() {
        withTempDirectory("work") { dir ->
            withTempFile("track", ".bed.gz", dir) { path ->
                withTempFile("control", ".bed.gz", dir) { control ->
                    // NOTE[oshpynov] we use .bed.gz here for the ease of sampling result save
                    SpanCLALongTest.sampleCoverage(path, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = true)
                    SpanCLALongTest.sampleCoverage(control, SpanCLALongTest.TO, SpanCLALongTest.BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    SpanCLA.main(arrayOf(
                        "analyze",
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

}