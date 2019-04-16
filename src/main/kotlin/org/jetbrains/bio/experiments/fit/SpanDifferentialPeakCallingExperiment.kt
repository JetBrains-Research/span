package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.span.Peak
import org.jetbrains.bio.span.getChromosomePeaks
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.MultiLabels
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.state.ZLHID
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths1: List<Pair<Path, Path?>>,
        paths2: List<Pair<Path, Path?>>,
        fragment: Fragment,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>
): SpanModelFitExperiment<Model, State>(
    genomeQuery,
    paths1 + paths2,
    MultiLabels.generate(TRACK1_PREFIX, paths1.size).toList() +
            MultiLabels.generate(TRACK2_PREFIX, paths2.size).toList(),
    fragment, binSize, modelFitter, modelClass, states, nullHypothesis
) {

    constructor(
            genomeQuery: GenomeQuery,
            paths1: Pair<Path, Path?>,
            paths2: Pair<Path, Path?>,
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
            reduceIds(paths1.flatMap { listOfNotNull(it.first, it.second) }.map { it.stemGz } +
                    listOf("vs") +
                    paths2.flatMap { listOfNotNull(it.first, it.second) }.map { it.stemGz }
                    + listOfNotNull(fragment.nullableInt, binSize).map { it.toString() })


    fun computeDirectedDifferencePeaks(fdr: Double,
                                       gap: Int): Pair<List<Peak>, List<Peak>> {
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            results.getChromosomePeaks(chromosome, fdr, gap, getData(chromosome))
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
                paths1: List<Pair<Path, Path?>>,
                paths2: List<Pair<Path, Path?>>,
                bin: Int,
                fragment: Fragment = AutoFragment
        ): SpanDifferentialPeakCallingExperiment<MLConstrainedNBHMM, ZLHID> {
            check(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            return if (paths1.size == 1 && paths2.size == 1) {
                SpanDifferentialPeakCallingExperiment(
                    genomeQuery, paths1.first(), paths2.first(),
                    fragment, bin,
                    semanticCheck(MLConstrainedNBHMM.fitter(1, 1), 1, 1),
                    MLConstrainedNBHMM::class.java,
                    ZLHID.values(), NullHypothesis.of(ZLHID.same())
                )
            } else {
                SpanDifferentialPeakCallingExperiment(
                    genomeQuery, paths1, paths2,
                    fragment, bin,
                    semanticCheck(MLConstrainedNBHMM.fitter(paths1.size, paths2.size), paths1.size, paths2.size),
                    MLConstrainedNBHMM::class.java,
                    ZLHID.values(), NullHypothesis.of(ZLHID.same())
                )
            }
        }

        private fun semanticCheck(fitter: Fitter<MLConstrainedNBHMM>, tracks1: Int, tracks2: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> by fitter {
                override fun fit(preprocessed: Preprocessed<DataFrame>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int,
                                 attempt: Int): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary(tracks1, tracks2)
                        }

                override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int,
                                 attempt: Int): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary(tracks1, tracks2)
                        }
            }
        }


    }
}