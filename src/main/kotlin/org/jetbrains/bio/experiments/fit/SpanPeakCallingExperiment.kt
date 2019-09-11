package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
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
) : SpanModelFitExperiment<Model, State>(
    genomeQuery,
    paths, null,
    MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
    fragment, binSize,
    modelFitter, modelClass, states, nullHypothesis, unique, fixedModelPath
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