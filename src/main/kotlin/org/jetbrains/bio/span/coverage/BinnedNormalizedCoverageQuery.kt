package org.jetbrains.bio.span.coverage

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.query.CachingQuery

/**
 * Binned version of [NormalizedCoverageQuery].
 */
class BinnedNormalizedCoverageQuery(
    private val normalizedCoverageQuery: NormalizedCoverageQuery,
    val binSize: Int
) : CachingQuery<Chromosome, IntArray>() {


    override val id: String
        get() = "binned_${binSize}_${normalizedCoverageQuery.id}"

    override val description: String
        get() = "Binned bin:${binSize}, coverage for ${normalizedCoverageQuery.description}"

    override fun getUncached(input: Chromosome): IntArray {
        return binnedCoverage[input]
    }

    private val binnedCoverage: GenomeMap<IntArray> by lazy {
        return@lazy genomeMap(normalizedCoverageQuery.genomeQuery, parallel = true) {
            it.range.slice(binSize).mapToInt { range ->
                normalizedCoverageQuery.apply(range.on(it))
            }.toArray()
        }
    }

}

fun List<BinnedNormalizedCoverageQuery>.binnedCoverageDataFrame(
    chromosome: Chromosome,
    labels: Array<String>
): DataFrame {
    var res = DataFrame()
    forEachIndexed { d, inputQuery ->
        val binnedCoverage = inputQuery.apply(chromosome)
        res = res.with(labels[d], binnedCoverage)
    }
    return res
}
