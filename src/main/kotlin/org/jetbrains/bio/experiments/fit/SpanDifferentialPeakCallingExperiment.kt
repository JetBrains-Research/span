package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.GenomeScoresQuery
import org.jetbrains.bio.coverage.scoresDataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.span.Peak
import org.jetbrains.bio.span.getChromosomePeaks
import org.jetbrains.bio.statistics.*
import org.jetbrains.bio.statistics.data.DataFrame
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.state.ZLHID
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment<Model : ClassificationModel, State : Any> : SpanModelFitExperiment<Model, State> {

    constructor(genomeQuery: GenomeQuery,
                paths1: Pair<Path, Path?>,
                paths2: Pair<Path, Path?>,
                binSize: Int, fragment: Int?,
                modelFitter: Fitter<Model>,
                modelClass: Class<Model>,
                states: Array<State>, nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery,
                    createDataQuery(genomeQuery,
                            listOf(paths1), listOf(paths2), fragment, binSize),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    constructor(genomeQuery: GenomeQuery,
                paths1: List<Pair<Path, Path?>>,
                paths2: List<Pair<Path, Path?>>,
                binSize: Int, fragment: Int?,
                modelFitter: Fitter<Model>, modelClass: Class<Model>,
                states: Array<State>, nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery,
                    createDataQuery(genomeQuery, paths1, paths2, fragment, binSize),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

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
                if (states.getAsFloat(it.startOffset / binSize, ZLHID.D.name) >
                        states.getAsFloat(it.startOffset / binSize, ZLHID.I.name)) {
                    highLow.add(it)
                } else {
                    lowHigh.add(it)
                }
            }
        }
        return highLow to lowHigh
    }


    override val id: String
        get() = "${dataQuery.id}_diff"

    companion object {
        const val TRACK1_PREFIX = "track1_"
        const val TRACK2_PREFIX = "track2_"

        internal fun createDataQuery(genomeQuery: GenomeQuery,
                                     paths1: List<Pair<Path, Path?>>,
                                     paths2: List<Pair<Path, Path?>>,
                                     fragment: Int?, binSize: Int): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                val scores1 = paths1.map { GenomeScoresQuery(genomeQuery, it.first, it.second, fragment, binSize) }
                val scores2 = paths2.map { GenomeScoresQuery(genomeQuery, it.first, it.second, fragment, binSize) }

                override fun getUncached(input: Chromosome): DataFrame {
                    return DataFrame.columnBind(scores1.scoresDataFrame(input, MultiLabels.generate(TRACK1_PREFIX, paths1.size)),
                            scores2.scoresDataFrame(input, MultiLabels.generate(TRACK2_PREFIX, paths2.size)))
                }

                override val id: String
                    get() = reduceIds(scores1.map { it.id } + listOf("vs") + scores2.map { it.id })

                override val description: String
                    get() = scores1.joinToString(";") { it.description } +
                            " vs " +
                            scores2.joinToString(";") { it.description }
            }
        }


        /**
         * Creates experiment for model-based comparison of binned coverage tracks for given queries.
         * Not restricted for single query and constrained for multiple queries.
         *
         * @return experiment [SpanDifferentialPeakCallingExperiment]
         */
        fun getExperiment(genomeQuery: GenomeQuery,
                          paths1: List<Pair<Path, Path?>>,
                          paths2: List<Pair<Path, Path?>>,
                          bin: Int,
                          fragment: Int? = null):
                SpanDifferentialPeakCallingExperiment<MLConstrainedNBHMM, ZLHID> {
            check(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            return if (paths1.size == 1 && paths2.size == 1) {
                SpanDifferentialPeakCallingExperiment(genomeQuery, paths1.first(), paths2.first(),
                        bin, fragment,
                        semanticCheck(MLConstrainedNBHMM.fitter(1, 1), 1, 1),
                        MLConstrainedNBHMM::class.java,
                        ZLHID.values(), NullHypothesis.of(ZLHID.same()))
            } else {
                SpanDifferentialPeakCallingExperiment(genomeQuery, paths1, paths2,
                        bin, fragment,
                        semanticCheck(MLConstrainedNBHMM.fitter(paths1.size, paths2.size), paths1.size, paths2.size),
                        MLConstrainedNBHMM::class.java,
                        ZLHID.values(), NullHypothesis.of(ZLHID.same()))
            }
        }

        private fun semanticCheck(fitter: Fitter<MLConstrainedNBHMM>, tracks1: Int, tracks2: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> by fitter {
                override fun fit(preprocessed: Preprocessed<DataFrame>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter).apply {
                            flipStatesIfNecessary(tracks1, tracks2)
                        }

                override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                                 title: String,
                                 threshold: Double,
                                 maxIter: Int): MLConstrainedNBHMM =
                        fitter.fit(preprocessed, title, threshold, maxIter).apply {
                            flipStatesIfNecessary(tracks1, tracks2)
                        }
            }
        }


    }
}