package org.jetbrains.bio.span

import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.coverage.SingleEndCoverage
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertTrue

class CoverageScoresQueryTest {
    private var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
    private var chromosome1: Chromosome = genomeQuery.get()[0]

    @Test
    fun testScoresWithControl() {
        val cond = SingleEndCoverage.builder(genomeQuery)
                .putAll(chromosome1, Strand.PLUS, 1, 2, 3, 4, 5, 10, 11, 15)
                .build(unique = false)

        val control = SingleEndCoverage.builder(genomeQuery)
                .putAll(chromosome1, Strand.PLUS, 0, 2, 4, 6, 10, 12, 14, 20, 21, 22, 25)
                .build(unique = false)

        assertTrue(Precision.equals(0.72, CoverageScoresQuery.computeScale(genomeQuery, cond, control), 0.01))

        assertEquals(listOf(0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1),
                CoverageScoresQuery.computeScores(chromosome1, cond, null, 1, 0.0).toList().subList(0, 16))

        assertEquals(listOf(0, 2, 0, 1, 0, 1, 0, 0, 0, 0),
                CoverageScoresQuery.computeScores(chromosome1, cond, control, 3, 1.0).toList().subList(0, 10))

        assertEquals(listOf(2, 0, 0, 1, 0),
                CoverageScoresQuery.computeScores(chromosome1, cond, control, 5, 0.5).toList().subList(0, 5))
    }

    @Test
    fun testScoresWithControlYD18H3K4me3LowQuality() {
        // This is the real data for low quality ChIP-Seq YD18 H3K4me3 Hg19 Chromosome1, first 100 bins x 200bp
        val treatmentPlus = intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 7, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        val treatmentMinus = intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 5, 3, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        val controlPlus = intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 1, 4, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
        val controlMinus = intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 3, 0, 1, 0, 0, 0, 0, 7, 2, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        val scale = 1.0

        // Scores from version 0.6.0.3800 with which best results were obtained.
        val expectedScores = intArrayOf(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

        val bin = 200


        val treatmentBuilder = SingleEndCoverage.builder(genomeQuery)
        treatmentPlus.forEachIndexed { index, v ->
            (0.until(v)).forEach { x ->
                treatmentBuilder.process(Location(index * bin + x, index * bin + x + 1, chromosome1, Strand.PLUS))
            }
        }
        treatmentMinus.forEachIndexed { index, v ->
            (0.until(v)).forEach { x ->
                treatmentBuilder.process(Location(index * bin + x, index * bin + x + 1, chromosome1, Strand.MINUS))
            }
        }
        val controlBuilder = SingleEndCoverage.builder(genomeQuery)
        controlPlus.forEachIndexed { index, v ->
            (0.until(v)).forEach { x ->
                controlBuilder.process(Location(index * bin + x, index * bin + x + 1, chromosome1, Strand.PLUS))
            }
        }
        controlMinus.forEachIndexed { index, v ->
            (0.until(v)).forEach { x ->
                controlBuilder.process(Location(index * bin + x, index * bin + x + 1, chromosome1, Strand.MINUS))
            }
        }

        assertEquals(expectedScores.toList(),
                CoverageScoresQuery.computeScores(chromosome1,
                        treatmentBuilder.build(unique = true),
                        controlBuilder.build(unique = true), bin, scale).toList().take(100))

    }

    @Test
    fun testControlMin0() {
        val cond = SingleEndCoverage.builder(genomeQuery)
                .putAll(chromosome1, Strand.PLUS, 1, 2, 4, 5, 10, 11, 15)
                .build(unique = false)

        val control = SingleEndCoverage.builder(genomeQuery)
                .putAll(chromosome1, Strand.PLUS, 1, 2, 3, 4, 5, 11, 15)
                .build(unique = false)
        assertEquals(listOf(0, 0, 1, 0),
                CoverageScoresQuery.computeScores(chromosome1, cond, control, 5, 1.0).toList().subList(0, 4))
    }
}

private fun SingleEndCoverage.Builder.putAll(chromosome: Chromosome,
                                    strand: Strand,
                                    vararg offsets: Int): SingleEndCoverage.Builder {
    data[chromosome, strand].addAll(offsets)
    return this
}
