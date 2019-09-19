package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.span.CoverageScoresQuery
import org.jetbrains.bio.span.Peak
import org.jetbrains.bio.span.getChromosomePeaks
import org.jetbrains.bio.span.scoresDataFrame
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.MultiLabels
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.state.ZLHID

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths1: List<SpanDataPaths>,
        paths2: List<SpanDataPaths>,
        fragment: Fragment,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>
): SpanModelFitExperiment<Model, Span1CompareFitInformation, State>(
    createEffectiveQueries(
        genomeQuery, paths1 + paths2,
        MultiLabels.generate(TRACK1_PREFIX, paths1.size).toList() +
                MultiLabels.generate(TRACK2_PREFIX, paths2.size).toList(),
        fragment, binSize
    ), Span1CompareFitInformation(
        genomeQuery,
        paths1, paths2,
        MultiLabels.generate(TRACK1_PREFIX, paths1.size).toList(),
        MultiLabels.generate(TRACK2_PREFIX, paths2.size).toList(),
        fragment, true, binSize
    ),
    modelFitter, modelClass, states, nullHypothesis
) {

    constructor(
            genomeQuery: GenomeQuery,
            paths1: SpanDataPaths,
            paths2: SpanDataPaths,
            fragment: Fragment, binSize: Int,
            modelFitter: Fitter<Model>,
            modelClass: Class<Model>,
            states: Array<State>, nullHypothesis: NullHypothesis<State>
    ): this(
        genomeQuery, listOf(paths1), listOf(paths2),
        fragment, binSize,
        modelFitter, modelClass, states, nullHypothesis
    )

    override val id: String =
            reduceIds(paths1.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOf("vs") +
                    paths2.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz }
                    + listOfNotNull(fragment.nullableInt, binSize).map { it.toString() })


    fun computeDirectedDifferencePeaks(fdr: Double,
            gap: Int): Pair<List<Peak>, List<Peak>> {
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            results.getChromosomePeaks(chromosome, fdr, gap, dataQuery.apply(chromosome))
        }
        val highLow = arrayListOf<Peak>()
        val lowHigh = arrayListOf<Peak>()
        genomeQuery.get().forEach { chromosome ->
            val states = getStatesDataFrame(chromosome)
            map[chromosome].forEach {
                if (states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.D.name) >
                        states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.I.name)) {
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
         * Not restricted for single query and constrained for multiple queries.
         *
         * @return experiment [SpanDifferentialPeakCallingExperiment]
         */
        fun getExperiment(
                genomeQuery: GenomeQuery,
                paths1: List<SpanDataPaths>,
                paths2: List<SpanDataPaths>,
                bin: Int,
                fragment: Fragment = AutoFragment
        ): SpanDifferentialPeakCallingExperiment<MLConstrainedNBHMM, ZLHID> {
            check(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            return if (paths1.size == 1 && paths2.size == 1) {
                SpanDifferentialPeakCallingExperiment(
                    genomeQuery, paths1.first(), paths2.first(),
                    fragment, bin,
                    MLConstrainedNBHMM.fitter(1, 1),
                    MLConstrainedNBHMM::class.java,
                    ZLHID.values(), NullHypothesis.of(ZLHID.same())
                )
            } else {
                SpanDifferentialPeakCallingExperiment(
                    genomeQuery, paths1, paths2,
                    fragment, bin,
                    MLConstrainedNBHMM.fitter(paths1.size, paths2.size),
                    MLConstrainedNBHMM::class.java,
                    ZLHID.values(), NullHypothesis.of(ZLHID.same())
                )
            }
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
) : FitInformation {

    constructor(
            genomeQuery: GenomeQuery,
            paths1: List<SpanDataPaths>,
            paths2: List<SpanDataPaths>,
            labels1: List<String>,
            labels2: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
    ): this(
        genomeQuery.build, paths1, paths2,
        labels1, labels2, fragment, unique, binSize,
        FitInformation.chromSizes(genomeQuery)
    )


    override fun scoresDataFrame(): Map<Chromosome, DataFrame> {
        val gq = genomeQuery()
        val queries1 = data1.map {
            CoverageScoresQuery(gq, it.treatment, it.control, fragment, binSize, unique)
        }
        val queries2 = data2.map {
            CoverageScoresQuery(gq, it.treatment, it.control, fragment, binSize, unique)
        }
        if (queries1.any { !it.ready } || queries2.any { !it.ready }) {
            return emptyMap()
        }
        return gq.get().associateBy({it}) {
            DataFrame.columnBind(
                queries1.scoresDataFrame(it, labels1.toTypedArray()),
                queries2.scoresDataFrame(it, labels2.toTypedArray())
            )
        }
    }

    companion object {
        const val VERSION: Int = 3
    }
}