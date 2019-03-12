package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.MultiLabels
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.hmm.MLAbstractHMM
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.state.ZLH
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths: List<Pair<Path, Path?>>,
        fragment: Int?,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>,
        unique: Boolean = true
) : SpanModelFitExperiment<Model, State>(
    genomeQuery,
    paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
    fragment, binSize,
    modelFitter, modelClass, states, nullHypothesis, unique
) {

    constructor(
            genomeQuery: GenomeQuery,
            paths: Pair<Path, Path?>,
            modelFitter: Fitter<Model>,
            modelClass: Class<Model>,
            fragment: Int?,
            binSize: Int,
            states: Array<State>,
            nullHypothesis: NullHypothesis<State>,
            unique: Boolean = true
    ): this(
        genomeQuery, listOf(paths),
        fragment, binSize,
        modelFitter, modelClass, states, nullHypothesis, unique
    )

    override val id: String =
            reduceIds(
                paths.flatMap { listOfNotNull(it.first, it.second) }.map { it.stemGz } +
                        listOfNotNull(fragment, binSize).map { it.toString() })


    companion object {
        const val TRACK_PREFIX = "track_"

        /**
         * Creates experiment for model-based enrichment of binned coverage tracks (e.g. ChIP-seq tracks)
         * for given number of [paths].
         * Not restricted for single query and constrained for multiple paths.
         *
         * @return experiment [SpanPeakCallingExperiment]
         */
        fun getExperiment(
                genomeQuery: GenomeQuery,
                paths: List<Pair<Path, Path?>>,
                bin: Int,
                fragment: Int? = null,
                unique: Boolean = true
        ): SpanPeakCallingExperiment<out MLAbstractHMM, ZLH> {
            check(paths.isNotEmpty()) { "No data" }
            return if (paths.size == 1) {
                SpanPeakCallingExperiment(
                    genomeQuery, paths.first(),
                    semanticCheck(MLFreeNBHMM.fitter()),
                    MLFreeNBHMM::class.java,
                    fragment, bin,
                    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                    unique
                )
            } else {
                SpanPeakCallingExperiment(
                    genomeQuery, paths,
                    fragment,
                    bin,
                    semanticCheck(MLConstrainedNBHMM.fitter(paths.size), paths.size),
                    MLConstrainedNBHMM::class.java,
                    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                    unique
                )
            }
        }

        private fun semanticCheck(fitter: Fitter<MLFreeNBHMM>): Fitter<MLFreeNBHMM> {
            return object : Fitter<MLFreeNBHMM> by fitter {

                override fun fit(preprocessed: Preprocessed<DataFrame>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int): MLFreeNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter).apply {
                            flipStatesIfNecessary()
                        }

                override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int): MLFreeNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter).apply {
                            flipStatesIfNecessary()
                        }
            }
        }


        private fun semanticCheck(fitter: Fitter<MLConstrainedNBHMM>, tracks: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> by fitter {

                override fun fit(preprocessed: Preprocessed<DataFrame>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter).apply {
                            flipStatesIfNecessary(tracks)
                        }

                override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter).apply {
                            flipStatesIfNecessary(tracks)
                        }
            }
        }


    }
}

