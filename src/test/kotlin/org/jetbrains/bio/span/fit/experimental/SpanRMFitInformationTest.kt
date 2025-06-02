package org.jetbrains.bio.span.fit.experimental

import com.google.common.math.IntMath
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.toPath
import org.jetbrains.bio.util.withTempFile
import org.junit.Assert.assertEquals
import org.junit.Test
import java.math.RoundingMode
import kotlin.test.assertTrue

class SpanRMFitInformationTest {

    val to1 = Genome["to1"].toQuery()

    @Test
    fun checkLoad() {
        val gq = GenomeQuery(Genome["to1"])
        val info = SpanRegrMixtureAnalyzeFitInformation(
            gq, SpanDataPaths("path_to_file".toPath(), "path_to_control".toPath()),
            null, "mapability.bigWig".toPath(), FixedFragment(42), false, 200
        )
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use {
                it.write(
                    """{
  "build": "to1",
  "paths": [
    {
      "treatment": "path_to_file",
      "control": "path_to_control"
    }
  ],
  "mapability_path": "mapability.bigWig",
  "labels": [
    "treatment_control"
  ],
  "fragment": 42,
  "unique": false,
  "bin_size": 200,
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrM": 10000,
    "chrX": 1000000
  },
  "fit.information.fqn": "org.jetbrains.bio.span.fit.experimental.SpanRegrMixtureAnalyzeFitInformation",  
  "version": 5
}"""
                )
            }
            assertEquals(info, SpanFitInformation.load(path))
        }
    }

    @Test
    fun testMapability() {
        val path = to1.genome.chromSizesPath.parent / "mapability.bigWig"
        val chrX = to1["chrX"]!!
        val binSize = 200
        /* for tests purposes we only generate mapability data for chrX */
        val mapabilityX = SpanRegrMixtureAnalyzeFitInformation.binnedMapability(chrX, path, binSize)
        sanityCheck(mapabilityX, chrX, binSize)
        val genomeMeanMapability = mapabilityX.sum() / mapabilityX.size
        assertTrue(genomeMeanMapability > 0.0, "Total mean mapability was 0.0")
        /*
            Genome mean mapability in the production code is calculated via BigWig summary method,
            which is not precise and depends on zoom levels.
        */
        val precision = 100.0 / chrX.length
        to1.get().filter { it != chrX }.forEach { chr ->
            val mapability = SpanRegrMixtureAnalyzeFitInformation.binnedMapability(chr, path, binSize)
            sanityCheck(mapability, chr, binSize)
            mapability.forEachIndexed { index, value ->
                assertEquals(
                    "In absence of data, mean mapability $value on $chr at bin $index should be equal " +
                            "to genome mean mapability $genomeMeanMapability",
                    genomeMeanMapability, value, precision
                )
            }
        }
    }

    companion object {
        fun sanityCheck(mapability: DoubleArray, chromosome: Chromosome, binSize: Int) {
            assertEquals(IntMath.divide(chromosome.length, binSize, RoundingMode.CEILING), mapability.size)
            mapability.forEachIndexed { index, value ->
                assertTrue(
                    value >= 0.0,
                    "Mean mapability value $value was negative for $chromosome at bin $index"
                )
                assertTrue(
                    value <= 1.0,
                    "Mean mapability value $value was greater than 1 for $chromosome at bin $index"
                )
            }
        }
    }

}