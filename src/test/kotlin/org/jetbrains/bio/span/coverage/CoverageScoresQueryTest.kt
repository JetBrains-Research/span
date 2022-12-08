package org.jetbrains.bio.span.coverage

import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.coverage.SingleEndCoverage
import org.junit.Test
import kotlin.test.assertTrue

class CoverageScoresQueryTest {
    private var genomeQuery: GenomeQuery = GenomeQuery(Genome["to1"])
    private var chromosome1: Chromosome = genomeQuery.get()[0]

    @Test
    fun testScoresWithControl() {
        val cond = SingleEndCoverage.builder(genomeQuery).apply {
            data[chromosome1, Strand.PLUS].addAll(intArrayOf(1, 2, 3, 4, 5, 10, 11, 15))
        }.build(unique = false)

        val control = SingleEndCoverage.builder(genomeQuery).apply {
            data[chromosome1, Strand.PLUS].addAll(intArrayOf(0, 2, 4, 6, 10, 12, 14, 20, 21, 22, 25))
        }.build(unique = false)

        val (scaleTreatment, scaleControl, beta) =
            NormalizedCoverageQuery.analyzeCoverage(genomeQuery, cond, control, 200)
        assertTrue(Precision.equals(1.375, scaleTreatment, 0.1))
        assertTrue(Precision.equals(1.0, scaleControl, 0.1))
        assertTrue(Precision.equals(0.63, beta, 0.01))
    }
}


