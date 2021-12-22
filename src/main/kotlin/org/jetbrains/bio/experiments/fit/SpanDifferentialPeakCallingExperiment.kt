package org.jetbrains.bio.experiments.fit

import com.google.common.math.DoubleMath.log2
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanFitInformation.Companion.chromSizes
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.coverage.BinnedCoverageScoresQuery
import org.jetbrains.bio.span.coverage.CoverageScoresQuery
import org.jetbrains.bio.span.coverage.binsScoresDataFrame
import org.jetbrains.bio.span.peaks.getChromosomeFdrGapPeaks
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.statistics.state.ZLHID
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import java.nio.file.Path

/**
 * Corresponds to Span `compare` invocation.
 *
 * The treatment-control pairs are split into two sets that are to compare.
 *
 * For each treatment-control pair, we compute binned DiffBind-like scores (see [BinnedCoverageScoresQuery] for details).
 * These scores are used as the input for a five-state multidimensional negative binomial HMM.
 * For each dimension `d`, there are two negative binomial distributions, low_d and high_d.
 * - ZERO state corresponds to zero emissions for all dimensions
 * - LOW state employs `low_d` emission for each dimension `d`
 * - HIGH state employs `high_d` emission for each dimension `d`
 * - INCREASED state employs `low_d` emission for each dimension `d` from the first set and `high_d` for the second set
 * - DECREASED state employs `high_d` emission for each dimension `d` from the first set and `low_d` for the second set
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment private constructor(
    fitInformation: Span1CompareFitInformation,
    threshold: Double,
    maxIter: Int
) : SpanModelFitExperiment<MLConstrainedNBHMM, Span1CompareFitInformation, ZLHID>(
    fitInformation,
    MLConstrainedNBHMM.fitter(fitInformation.data1.size, fitInformation.data2.size),
    MLConstrainedNBHMM::class.java,
    ZLHID.values(), NullHypothesis.of(ZLHID.same()),
    threshold = threshold,
    maxIter = maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

    fun computeDirectedDifferencePeaks(
        fdr: Double,
        gap: Int
    ): Pair<List<Peak>, List<Peak>> {
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            results.getChromosomeFdrGapPeaks(chromosome, fdr, gap)
        }
        val highLow = arrayListOf<Peak>()
        val lowHigh = arrayListOf<Peak>()
        genomeQuery.get().forEach { chromosome ->
            val states = getStatesDataFrame(chromosome)
            map[chromosome].forEach {
                if (states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.D.name) >
                    states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.I.name)
                ) {
                    highLow.add(it)
                } else {
                    lowHigh.add(it)
                }
            }
        }
        return highLow to lowHigh
    }


    companion object {
        const val TRACK1_PREFIX = "track1_"
        const val TRACK2_PREFIX = "track2_"

        /**
         * Creates experiment for model-based comparison of binned coverage tracks for given queries.
         *
         * @return experiment [SpanDifferentialPeakCallingExperiment]
         */
        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths1: List<SpanDataPaths>,
            paths2: List<SpanDataPaths>,
            bin: Int,
            fragment: Fragment,
            unique: Boolean,
            threshold: Double,
            maxIter: Int
        ): SpanDifferentialPeakCallingExperiment {
            check(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            val fitInformation = Span1CompareFitInformation.effective(
                genomeQuery,
                paths1, paths2,
                MultiLabels.generate(TRACK1_PREFIX, paths1.size).toList(),
                MultiLabels.generate(TRACK2_PREFIX, paths2.size).toList(),
                fragment, unique, bin
            )
            return SpanDifferentialPeakCallingExperiment(fitInformation, threshold, maxIter)
        }
    }
}

data class Span1CompareFitInformation(
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
            val data = data1 + data2
            val labels = labels1 + labels2
            return object : CachingQuery<Chromosome, DataFrame>() {
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
        }

    @Transient
    private lateinit var scoreQueries1: List<CoverageScoresQuery>

    @Transient
    private lateinit var scoreQueries2: List<CoverageScoresQuery>

    override fun prepareScores() {
        scoreQueries1 = data1.map {
            CoverageScoresQuery(genomeQuery(), it.treatment, it.control, fragment, unique, showLibraryInfo = false)
        }
        scoreQueries2 = data2.map {
            CoverageScoresQuery(genomeQuery(), it.treatment, it.control, fragment, unique, showLibraryInfo = false)
        }
    }

    /**
     * Return log2 fold change of average summary coverage across data
     */
    override fun score(chromosomeRange: ChromosomeRange): Double {
        return if (scoreQueries1.all { it.ready } && scoreQueries2.all { it.ready }) {
            val score1 = scoreQueries1.sumOf { it.apply(chromosomeRange) }.toDouble() / scoreQueries1.size
            val score2 = scoreQueries2.sumOf { it.apply(chromosomeRange) }.toDouble() / scoreQueries2.size
            if (score2 != 0.0) log2(score1) - log2(score2) else Double.MAX_VALUE
        } else {
            0.0
        }
    }

    companion object {
        const val VERSION: Int = 3

        fun effective(
            genomeQuery: GenomeQuery,
            paths1: List<SpanDataPaths>,
            paths2: List<SpanDataPaths>,
            labels1: List<String>,
            labels2: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
        ): Span1CompareFitInformation {
            return Span1CompareFitInformation(
                genomeQuery.build, paths1, paths2, labels1, labels2, fragment, unique, binSize,
                chromSizes(
                    SpanModelFitExperiment.effectiveGenomeQuery(
                        genomeQuery, paths1 + paths2, fragment, unique
                    )
                )
            )
        }
    }
}