package org.jetbrains.bio.span.fit

import com.google.common.math.DoubleMath
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

data class SpanCompareFitInformation(
    override val build: String,
    val data1: List<SpanDataPaths>,
    val data2: List<SpanDataPaths>,
    val labels1: List<String>,
    val labels2: List<String>,
    val fragment: Fragment,
    val unique: Boolean,
    override val binSize: Int,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : SpanFitInformation {

    override val id
        get() = reduceIds(
            data1.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOf("vs") +
                    data2.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
        )

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            prepareData()
            val binnedScoresQueries = (scoreQueries1!! + scoreQueries2!!).map { BinnedCoverageScoresQuery(it, binSize) }
            val labels = labels1 + labels2
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return binnedScoresQueries.binsScoresDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(binnedScoresQueries.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
            }
        }

    @Transient
    var scoreQueries1: List<CoverageScoresQuery>? = null

    @Transient
    var scoreQueries2: List<CoverageScoresQuery>? = null

    @Synchronized
    override fun prepareData() {
        if (scoreQueries1 == null) {
            scoreQueries1 = data1.map {
                CoverageScoresQuery(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    fragment,
                    unique,
                    showLibraryInfo = false
                )
            }
        }
        if (scoreQueries2 == null) {
            scoreQueries2 = data2.map {
                CoverageScoresQuery(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    fragment,
                    unique,
                    showLibraryInfo = false
                )
            }
        }
    }

    /**
     * Return log2 fold change of average summary coverage across data
     */
    @Synchronized
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(scoreQueries1 != null && scoreQueries2 != null) {
            "Please use prepareData before!"
        }

        return if (scoreQueries1!!.all { it.ready } && scoreQueries2!!.all { it.ready }) {
            val score1 = scoreQueries1!!.sumOf { it.apply(chromosomeRange) }.toDouble() / scoreQueries1!!.size
            val score2 = scoreQueries2!!.sumOf { it.apply(chromosomeRange) }.toDouble() / scoreQueries2!!.size
            if (score2 != 0.0) DoubleMath.log2(score1) - DoubleMath.log2(score2) else Double.MAX_VALUE
        } else {
            0.0
        }
    }

    @Synchronized
    override fun scaledTreatmentScore(chromosomeRange: ChromosomeRange): Double {
        check(scoreQueries1 != null && scoreQueries2 != null) {
            "Please use prepareData before!"
        }
        return scoreQueries1!!.sumOf { it.scaledTreatment(chromosomeRange) } / scoreQueries1!!.size
    }

    @Synchronized
    override fun scaledControlScore(chromosomeRange: ChromosomeRange): Double {
        check(scoreQueries1 != null && scoreQueries2 != null) {
            "Please use prepareData before!"
        }
        return scoreQueries2!!.sumOf { it.scaledTreatment(chromosomeRange) } / scoreQueries2!!.size
    }

    override fun hasControlData(): Boolean {
        check(scoreQueries1 != null && scoreQueries2 != null) {
            "Please use prepareData before!"
        }
        return true
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 4

        fun effective(
            genomeQuery: GenomeQuery,
            paths1: List<SpanDataPaths>,
            paths2: List<SpanDataPaths>,
            labels1: List<String>,
            labels2: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
        ): SpanCompareFitInformation {
            return SpanCompareFitInformation(
                genomeQuery.build, paths1, paths2, labels1, labels2, fragment, unique, binSize,
                SpanFitInformation.chromSizes(
                    SpanModelFitExperiment.filterGenomeQueryWithData(
                        genomeQuery, paths1 + paths2, fragment, unique
                    )
                )
            )
        }
    }
}