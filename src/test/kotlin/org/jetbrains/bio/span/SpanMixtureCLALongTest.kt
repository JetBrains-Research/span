package org.jetbrains.bio.span

import kotlinx.support.jdk7.use
import org.jetbrains.bio.Tests.assertIn
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.span.SpanCLALongTest.Companion.THREADS
import org.jetbrains.bio.span.SpanCLALongTest.Companion.TO
import org.jetbrains.bio.span.coverage.SpanCoverageSampler.sampleCoverage
import org.jetbrains.bio.span.fit.SpanModelType
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.*
import org.junit.After
import org.junit.Before
import org.junit.Test
import java.io.FileReader
import java.util.*
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class SpanMixtureCLALongTest {

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
                    sampleCoverage(path, TO, SpanPeakCallingExperiment.SPAN_DEFAULT_BIN, goodQuality = true)
                    sampleCoverage(control, TO, SpanPeakCallingExperiment.SPAN_DEFAULT_BIN, goodQuality = false)

                    val chromsizes = Genome["to1"].chromSizesPath.toString()
                    SpanCLA.main(
                        arrayOf(
                            "analyze",
                            "-cs", chromsizes,
                            "--workdir", dir.toString(),
                            "-t", path.toString(),
                            "-c", control.toString(),
                            "--threads", THREADS.toString(),
                            "--model-type", SpanModelType.NB2Z_MIXTURE.id
                        )
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
                                            "50"
                                        )
                                    )
                                }.${SpanModelType.NB2Z_MIXTURE.extension}"
                            ).size
                    )
                }
            }
        }
    }


    @Test
    fun analyzeNB2ZMixtureSampledEnrichment() {
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
                SpanPeakCallingExperiment.SPAN_DEFAULT_BIN,
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
                        "--debug",
                        "-cs", Genome["to1"].chromSizesPath.toString(),
                        "-w", dir.toString(),
                        "--peaks", bedPath.toString(),
                        "-fdr", "0.05",
                        "-t", path.toString(),
                        "--model-type", SpanModelType.NB2Z_MIXTURE.id,
                    )
                )
                // Check created bed file
                assertTrue(
                    Location(
                        1100 * SpanPeakCallingExperiment.SPAN_DEFAULT_BIN,
                        1900 * SpanPeakCallingExperiment.SPAN_DEFAULT_BIN,
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

}
