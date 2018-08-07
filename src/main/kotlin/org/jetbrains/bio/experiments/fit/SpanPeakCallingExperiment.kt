package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.GenomeScoresQuery
import org.jetbrains.bio.coverage.scoresDataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.statistics.*
import org.jetbrains.bio.statistics.data.DataFrame
import org.jetbrains.bio.statistics.hmm.MLAbstractHMM
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.state.ZLH
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel, State : Any> : SpanFitExperiment<Model, State> {

    constructor(genomeQuery: GenomeQuery,
                paths: Pair<Path, Path?>,
                modelFitter: Fitter<Model>,
                modelClass: Class<Model>,
                binSize: Int,
                fragment: Int?,
                states: Array<State>,
                nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery,
                    createDataQuery(genomeQuery, listOf(paths), binSize, fragment, arrayOf(X_PREFIX)),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    constructor(genomeQuery: GenomeQuery,
                paths: List<Pair<Path, Path?>>,
                modelFitter: Fitter<Model>,
                modelClass: Class<Model>,
                binSize: Int,
                fragment: Int?,
                states: Array<State>,
                nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery,
                    createDataQuery(genomeQuery, paths, binSize, fragment, MultiLabels.generate(D_PREFIX, paths.size)),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    override val id: String
        get() = dataQuery.id

    companion object {
        const val X_PREFIX = "x"
        const val D_PREFIX = "d"

        private fun createDataQuery(genomeQuery: GenomeQuery,
                                    paths: List<Pair<Path, Path?>>,
                                    binSize: Int,
                                    fragment: Int?,
                                    labels: Array<String>): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                val scores = paths.map { GenomeScoresQuery.create(genomeQuery, it.first, it.second, fragment, binSize) }

                override fun getUncached(input: Chromosome): DataFrame {
                    return scores.scoresDataFrame(input, labels)
                }

                override val id: String
                    get() = reduceIds(scores.map { it.id })

                override val description: String
                    get() = scores.joinToString(";") { it.description }
            }
        }

        /**
         * Creates experiment for model-based enrichment of binned coverage tracks (e.g. ChIP-seq tracks)
         * for given number of [paths].
         * Not restricted for single query and constrained for multiple paths.
         *
         * @return experiment [SpanPeakCallingExperiment]
         */
        fun getExperiment(genomeQuery: GenomeQuery,
                          paths: List<Pair<Path, Path?>>,
                          bin: Int,
                          fragment: Int? = null):
                SpanPeakCallingExperiment<out MLAbstractHMM, ZLH> {
            check(paths.isNotEmpty()) { "No data" }
            return if (paths.size == 1) {
                SpanPeakCallingExperiment(genomeQuery, paths.first(),
                        semanticCheck(MLFreeNBHMM.fitter()),
                        MLFreeNBHMM::class.java,
                        bin, fragment,
                        ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L))
            } else {
                SpanPeakCallingExperiment(genomeQuery, paths,
                        semanticCheck(MLConstrainedNBHMM.fitter(paths.size), paths.size),
                        MLConstrainedNBHMM::class.java,
                        bin, fragment,
                        ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L))
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

