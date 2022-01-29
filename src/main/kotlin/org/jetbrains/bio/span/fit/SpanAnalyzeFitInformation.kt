package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.span.coverage.BinnedCoverageScoresQuery
import org.jetbrains.bio.span.coverage.CoverageScoresQuery
import org.jetbrains.bio.span.coverage.binsScoresDataFrame
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz

/**
 * Since all the chromosomes are squashed in [SpanModelFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 *
 * [labels] refer to the coverage dataframe column labels, not to the supervised learning annotations.
 */
data class SpanAnalyzeFitInformation(
    override val build: String,
    override val data: List<SpanDataPaths>,
    val labels: List<String>,
    override val fragment: Fragment,
    override val unique: Boolean,
    override val binSize: Int,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : AbstractSpanAnalyzeFitInformation {

    constructor(
        genomeQuery: GenomeQuery,
        paths: List<SpanDataPaths>,
        labels: List<String>,
        fragment: Fragment,
        unique: Boolean,
        binSize: Int
    ) : this(
        genomeQuery.build, paths,
        labels, fragment, unique, binSize,
        SpanFitInformation.chromSizes(genomeQuery)
    )

    override val id: String
        get() = reduceIds(
            data.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
        )

    override val dataQuery: Query<Chromosome, DataFrame>
        get() = object : CachingQuery<Chromosome, DataFrame>() {
            val scores = data.map {
                BinnedCoverageScoresQuery(
                    CoverageScoresQuery(
                        genomeQuery(),
                        it.treatment,
                        it.control,
                        fragment,
                        unique
                    ), binSize
                )
            }

            override fun getUncached(input: Chromosome): DataFrame {
                return scores.binsScoresDataFrame(input, labels.toTypedArray())
            }

            override val id: String
                get() = reduceIds(scores.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
        }

    @Transient
    private var scoreQueries: List<CoverageScoresQuery>? = null

    override fun prepareScores() {
        scoreQueries = data.map {
            CoverageScoresQuery(genomeQuery(), it.treatment, it.control, fragment, unique, showLibraryInfo = false)
        }
    }

    /**
     * Returns summary coverage averaged by tracks
     */
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(scoreQueries != null) {
            "Please use prepareScores before!"
        }
        return if (scoreQueries!!.all { it.ready }) {
            scoreQueries!!.sumOf { it.apply(chromosomeRange) }.toDouble() / scoreQueries!!.size
        } else {
            0.0
        }
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 4

        fun createFitInformation(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            labels: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
        ): SpanAnalyzeFitInformation {
            val genomeQueryWithData =
                SpanModelFitExperiment.filterGenomeQueryWithData(genomeQuery, paths, fragment, unique)
            return SpanAnalyzeFitInformation(
                genomeQueryWithData.build,
                paths,
                labels,
                fragment,
                unique,
                binSize,
                SpanFitInformation.chromSizes(genomeQueryWithData)
            )
        }
    }
}
