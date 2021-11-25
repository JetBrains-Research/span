package org.jetbrains.bio.span.coverage

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.query.CachingQuery

/**
 * Coverage query to reproduce DiffBind-like scores for binned genome.
 *
 * See [CoverageScoresQuery] for realization details.
 */
class BinnedCoverageScoresQuery(
    val coverageScoresQuery: CoverageScoresQuery,
    val binSize: Int
) : CachingQuery<Chromosome, IntArray>() {


    override val id: String
        get() = "binned_${binSize}_${coverageScoresQuery.id}"

    override val description: String
        get() = "Binned bin:${binSize}, coverage for ${coverageScoresQuery.description}"

    override fun getUncached(input: Chromosome): IntArray {
        return scores[input]
    }

    val scores: GenomeMap<IntArray> by lazy {
        return@lazy genomeMap(coverageScoresQuery.genomeQuery, parallel = true) {
            computeBinnedScores(it, binSize)
        }
    }

    val ready: Boolean get() = coverageScoresQuery.ready

    /**
     * Returns the scores of a given [chromosome] sliced into binSizes with [binSize] width.
     * Score is precisely [Coverage] within bins if no control given, or DiffBind-like score.
     */
    internal fun computeBinnedScores(chromosome: Chromosome, binSize: Int): IntArray {
        return chromosome.range.slice(binSize).mapToInt { bin ->
            coverageScoresQuery.apply(bin.on(chromosome))
        }.toArray()
    }

}

fun List<BinnedCoverageScoresQuery>.binsScoresDataFrame(chromosome: Chromosome, labels: Array<String>): DataFrame {
    var res = DataFrame()
    forEachIndexed { d, inputQuery ->
        val binnedCoverage = inputQuery.apply(chromosome)
        res = res.with(labels[d], binnedCoverage)
    }
    return res
}
