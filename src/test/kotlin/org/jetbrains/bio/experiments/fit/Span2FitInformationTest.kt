package org.jetbrains.bio.experiments.fit

import com.google.common.math.IntMath
import org.jetbrains.bio.Tests
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.toQuery
import org.jetbrains.bio.util.div
import org.junit.Test
import java.math.RoundingMode
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class Span2FitInformationTest {

    val to1 = Genome["to1"].toQuery()

    @Test
    fun testMapability() {
        val path = to1.genome.chromSizesPath.parent / "mapability.bigWig"
        val chrX = to1["chrX"]!!
        val binSize = 200
        /* for tests purposes we only generate mapability data for chrX */
        val mapabilityX = Span2FitInformation.binnedMapability(chrX, path, binSize)
        sanityCheck(mapabilityX, chrX, binSize)
        val genomeMeanMapability = mapabilityX.sum() / mapabilityX.size
        assertTrue(genomeMeanMapability > 0.0, "Total mean mapability was 0.0")
        /**
         * Genome mean mapability in the production code is calculated via BigWig summary method,
         * which is not precise and depends on zoom levels.
         */
        val precision = 10.0 / chrX.length
        to1.get().filter { it != chrX }.forEach { chr ->
            val mapability = Span2FitInformation.binnedMapability(chr, path, binSize)
            sanityCheck(mapability, chr, binSize)
            mapability.forEachIndexed { index, value ->
                Tests.assertEquals(
                    genomeMeanMapability, value, precision,
                    "In absence of data, mean mapability $value on $chr at bin $index should be equal " +
                            "to genome mean mapability $genomeMeanMapability"
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