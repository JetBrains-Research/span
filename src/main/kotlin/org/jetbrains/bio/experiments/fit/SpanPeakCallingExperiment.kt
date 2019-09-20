package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanFitInformation.Companion.chromSizes
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.span.CoverageScoresQuery
import org.jetbrains.bio.span.scoresDataFrame
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.MultiLabels
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel, State : Any> private constructor(
        fitInformation: Span1AnalyzeFitInformation,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>,
        fixedModelPath: Path?
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, State>(
    fitInformation, modelFitter, modelClass, states, nullHypothesis, fixedModelPath
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

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
                paths: List<SpanDataPaths>,
                bin: Int,
                fragment: Fragment = AutoFragment,
                unique: Boolean = true,
                fixedModelPath: Path? = null
        ): SpanPeakCallingExperiment<out ClassificationModel, ZLH> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            return if (paths.size == 1) {
                SpanPeakCallingExperiment(
                    fitInformation,
                    MLFreeNBHMM.fitter().multiStarted(),
                    MLFreeNBHMM::class.java, ZLH.values(),
                    NullHypothesis.of(ZLH.Z, ZLH.L),
                    fixedModelPath
                )
            } else {
                SpanPeakCallingExperiment(
                    fitInformation,
                    MLConstrainedNBHMM.fitter(paths.size).multiStarted(),
                    MLConstrainedNBHMM::class.java, ZLH.values(),
                    NullHypothesis.of(ZLH.Z, ZLH.L),
                    fixedModelPath
                )
            }
        }
    }
}

/**
 * Since all the chromosomes are squashed in [SpanModelFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 *
 * [labels] refer to the coverage dataframe column labels, not to the supervised learning annotations.
 */
data class Span1AnalyzeFitInformation(
        override val build: String,
        override val data: List<SpanDataPaths>,
        override val labels: List<String>,
        override val fragment: Fragment,
        override val unique: Boolean,
        override val binSize: Int,
        override val chromosomesSizes: LinkedHashMap<String, Int>
) : SpanAnalyzeFitInformation {

    constructor(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            labels: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
    ): this(
        genomeQuery.build, paths,
        labels, fragment, unique, binSize,
        chromSizes(genomeQuery)
    )

    override val id: String =
            reduceIds(
                data.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                        listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
            )

    override val dataQuery: Query<Chromosome, DataFrame>
        get() = object : CachingQuery<Chromosome, DataFrame>() {
            val scores = data.map {
                CoverageScoresQuery(genomeQuery(), it.treatment, it.control, fragment, binSize, unique)
            }

            override fun getUncached(input: Chromosome): DataFrame {
                return scores.scoresDataFrame(input, labels.toTypedArray())
            }

            override val id: String
                get() = reduceIds(scores.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
        }


    override fun scoresDataFrame(): Map<Chromosome, DataFrame> {
        val gq = genomeQuery()
        val queries = data.map {
            CoverageScoresQuery(gq, it.treatment, it.control, fragment, binSize, unique)
        }
        if (queries.any { !it.ready }) {
            return emptyMap()
        }
        return gq.get().associateBy({it}) {
            queries.scoresDataFrame(it, labels.toTypedArray())
        }
    }

    companion object {
        const val VERSION: Int = 3

        fun effective(
                genomeQuery: GenomeQuery,
                paths: List<SpanDataPaths>,
                labels: List<String>,
                fragment: Fragment,
                unique: Boolean,
                binSize: Int
        ): Span1AnalyzeFitInformation {
            val effectiveGQ = SpanModelFitExperiment.effectiveGenomeQuery(genomeQuery, paths, fragment, unique)
            return Span1AnalyzeFitInformation(
                effectiveGQ.build, paths, labels, fragment, unique, binSize, chromSizes(effectiveGQ)
            )
        }
    }
}