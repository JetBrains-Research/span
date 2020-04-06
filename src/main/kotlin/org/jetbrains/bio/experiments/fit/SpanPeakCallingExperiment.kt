package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanFitInformation.Companion.chromSizes
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.stemGz
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
import org.jetbrains.bio.util.reduceIds
import java.nio.file.Path

/**
 * Corresponds to Span `analyze --type nbhmm` invocation.
 *
 * For each treatment-control pair, we compute binned DiffBind-like scores (see [CoverageScoresQuery] for details).
 * These scores are used as the input for a three-state multidimensional negative binomial HMM.
 * For each dimension `d`, there are two negative binomial distributions, low_d and high_d.
 * - ZERO state corresponds to zero emissions for all dimensions
 * - LOW state employs `low_d` emission for each dimension `d`
 * - HIGH state employs `high_d` emission for each dimension `d`
 *
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel> private constructor(
        fitInformation: Span1AnalyzeFitInformation,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        fixedModelPath: Path?,
        threshold: Double,
        maxIter: Int
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, ZLH>(
        fitInformation, modelFitter, modelClass, ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L), fixedModelPath,
        threshold, maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

    companion object {

        const val SPAN_DEFAULT_BIN = 200
        const val SPAN_DEFAULT_FDR = 1E-6
        const val SPAN_DEFAULT_GAP = 5

        const val SPAN_REPLICATED_DEFAULT_FDR = 1E-10

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
                fixedModelPath: Path? = null,
                threshold: Double = Fitter.THRESHOLD,
                maxIter: Int = Fitter.MAX_ITERATIONS,
                multistarts: Int = Fitter.MULTISTARTS,
                multistartIter: Int = Fitter.MULTISTART_ITERATIONS
        ): SpanPeakCallingExperiment<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                    genomeQuery, paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
                    fragment, unique, bin
            )
            return if (paths.size == 1) {
                SpanPeakCallingExperiment(
                        fitInformation,
                        when {
                            multistarts > 0 ->
                                MLFreeNBHMM.fitter().multiStarted(multistarts, multistartIter)
                            else ->
                                MLFreeNBHMM.fitter()
                        },
                        MLFreeNBHMM::class.java,
                        fixedModelPath,
                        threshold,
                        maxIter
                )
            } else {
                SpanPeakCallingExperiment(
                        fitInformation,
                        when {
                            multistarts > 0 ->
                                MLConstrainedNBHMM.fitter(paths.size).multiStarted(multistarts, multistartIter)
                            else ->
                                MLConstrainedNBHMM.fitter(paths.size)
                        },
                        MLConstrainedNBHMM::class.java,
                        fixedModelPath,
                        threshold,
                        maxIter
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
        val labels: List<String>,
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
    ) : this(
            genomeQuery.build, paths,
            labels, fragment, unique, binSize,
            chromSizes(genomeQuery)
    )

    override val id: String
        get() = reduceIds(
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
        return gq.get().associateBy({ it }) {
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
