package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.FitInformation.Companion.chromSizes
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
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
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths: List<SpanDataPaths>,
        fragment: Fragment,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>,
        unique: Boolean = true,
        fixedModelPath: Path? = null
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, State>(
    createEffectiveQueries(
        genomeQuery, paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(), fragment, binSize, unique
    ),
    Span1AnalyzeFitInformation(
        genomeQuery, paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
        fragment, unique, binSize
    ),
    modelFitter, modelClass, states, nullHypothesis, fixedModelPath
) {

    constructor(
            genomeQuery: GenomeQuery,
            paths: SpanDataPaths,
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
                paths.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
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
                paths: List<SpanDataPaths>,
                bin: Int,
                fragment: Fragment = AutoFragment,
                unique: Boolean = true,
                fixedModelPath: Path? = null
        ): SpanPeakCallingExperiment<out ClassificationModel, ZLH> {
            check(paths.isNotEmpty()) { "No data" }
            return if (paths.size == 1) {
                SpanPeakCallingExperiment(
                    genomeQuery, paths.first(),
                    MLFreeNBHMM.fitter().multiStarted(),
                    MLFreeNBHMM::class.java,
                    fragment, bin,
                    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                    unique,
                    fixedModelPath
                )
            } else {
                SpanPeakCallingExperiment(
                    genomeQuery, paths,
                    fragment,
                    bin,
                    MLConstrainedNBHMM.fitter(paths.size).multiStarted(),
                    MLConstrainedNBHMM::class.java,
                    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                    unique,
                    fixedModelPath
                )
            }
        }
    }
}

enum class SpanModel(val description: String) {
    NB_HMM("negative binomial HMM"), POISSON_REGRESSION_MIXTURE("Poisson regression mixture");

    override fun toString() = description
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
    }
}