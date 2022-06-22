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
        get() {
            prepareData()
            val binnedScoresQueries = scoreQueries!!.map { BinnedCoverageScoresQuery(it, binSize) }
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return binnedScoresQueries.binsScoresDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(binnedScoresQueries.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
            }
        }

    @Transient
    var scoreQueries: List<CoverageScoresQuery>? = null

    @Synchronized
    override fun prepareData() {
        if (scoreQueries == null) {
            scoreQueries = data.map {
                CoverageScoresQuery(genomeQuery(), it.treatment, it.control, fragment, unique, showLibraryInfo = false)
            }
        }
    }

    /**
     * Returns average coverage by tracks
     */
    @Synchronized
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(scoreQueries != null) {
            "Please use prepareData before!"
        }
        return if (scoreQueries!!.all { it.ready }) {
            scoreQueries!!.sumOf { it.apply(chromosomeRange) }.toDouble() / scoreQueries!!.size
        } else {
            0.0
        }
    }

    @Synchronized
    override fun scaledTreatmentScore(chromosomeRange: ChromosomeRange): Double {
        check(scoreQueries != null) {
            "Please use prepareData before!"
        }
        return scoreQueries!!.sumOf { it.scaledTreatment(chromosomeRange) } / scoreQueries!!.size
    }

    @Synchronized
    override fun scaledControlScore(chromosomeRange: ChromosomeRange): Double? {
        check(scoreQueries != null) {
            "Please use prepareData before!"
        }
        if (scoreQueries!!.any { it.controlReads == null }) {
            return null
        }
        return scoreQueries!!.sumOf { it.scaledControl(chromosomeRange)!! } / scoreQueries!!.size
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
