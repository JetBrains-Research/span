package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.MultiLabels
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.mixture.ZeroPoissonMixture
import org.jetbrains.bio.statistics.state.ZLH
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths: List<SpanPathsToData>,
        fragment: Fragment,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>,
        unique: Boolean = true,
        fixedModelPath: Path? = null
) : SpanModelFitExperiment<Model, State>(
    genomeQuery,
    paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
    fragment, binSize,
    modelFitter, modelClass, states, nullHypothesis, unique, fixedModelPath
) {

    constructor(
            genomeQuery: GenomeQuery,
            paths: SpanPathsToData,
            modelFitter: Fitter<Model>,
            modelClass: Class<Model>,
            fragment: Fragment,
            binSize: Int,
            states: Array<State>,
            nullHypothesis: NullHypothesis<State>,
            unique: Boolean = true,
            fixedModelPath: Path? = null
    ) : this(
        genomeQuery, listOf(paths),
        fragment, binSize,
        modelFitter, modelClass, states, nullHypothesis, unique, fixedModelPath
    )

    override val id: String =
            reduceIds(
                paths.flatMap { listOfNotNull(it.pathTreatment, it.pathInput) }.map { it.stemGz } +
                        listOfNotNull(fragment.nullableInt, binSize).map { it.toString() })


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
                paths: List<SpanPathsToData>,
                bin: Int,
                fragment: Fragment = AutoFragment,
                unique: Boolean = true,
                fixedModelPath: Path? = null,
                modelType: SpanModel = SpanModel.NB_HMM
        ): SpanPeakCallingExperiment<out ClassificationModel, ZLH> {
            check(paths.isNotEmpty()) { "No data" }
            return if (paths.size == 1) {
                if (modelType == SpanModel.NB_HMM) {
                    SpanPeakCallingExperiment(
                        genomeQuery, paths.first(),
                        semanticCheck(MLFreeNBHMM.fitter()),
                        MLFreeNBHMM::class.java,
                        fragment, bin,
                        ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                        unique,
                        fixedModelPath
                    )
                } else {
                    SpanPeakCallingExperiment(
                        genomeQuery, paths.first(),
                        semanticCheckZPM(ZeroPoissonMixture.fitter()),
                        ZeroPoissonMixture::class.java,
                        fragment, bin,
                        ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                        unique,
                        fixedModelPath
                    )
                }
            } else {
                SpanPeakCallingExperiment(
                    genomeQuery, paths,
                    fragment,
                    bin,
                    semanticCheck(MLConstrainedNBHMM.fitter(paths.size), paths.size).multiStarted(),
                    MLConstrainedNBHMM::class.java,
                    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                    unique,
                    fixedModelPath
                )
            }
        }

        private fun semanticCheckZPM(fitter: Fitter<ZeroPoissonMixture>): Fitter<ZeroPoissonMixture> {
            return object : Fitter<ZeroPoissonMixture> by fitter {

                override fun fit(
                        preprocessed: Preprocessed<DataFrame>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): ZeroPoissonMixture =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary()
                        }

                override fun fit(
                        preprocessed: List<Preprocessed<DataFrame>>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): ZeroPoissonMixture =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary()
                        }
            }
        }

        private fun semanticCheck(fitter: Fitter<MLFreeNBHMM>): Fitter<MLFreeNBHMM> {
            return object : Fitter<MLFreeNBHMM> by fitter {

                override fun fit(
                        preprocessed: Preprocessed<DataFrame>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): MLFreeNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary()
                        }

                override fun fit(
                        preprocessed: List<Preprocessed<DataFrame>>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): MLFreeNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary()
                        }
            }
        }


        private fun semanticCheck(fitter: Fitter<MLConstrainedNBHMM>, tracks: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> by fitter {

                override fun fit(
                        preprocessed: Preprocessed<DataFrame>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary(tracks)
                        }

                override fun fit(
                        preprocessed: List<Preprocessed<DataFrame>>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary(tracks)
                        }
            }
        }
    }
}

enum class SpanModel {
    NB_HMM, POISSON_REGRESSION_MIXTURE
}