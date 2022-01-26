package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.coverage.BinnedCoverageScoresQuery
import org.jetbrains.bio.span.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.span.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Span `analyze --type nbhmm` invocation.
 *
 * For each treatment-control pair, we compute binned DiffBind-like scores (see [BinnedCoverageScoresQuery] for details).
 * These scores are used as the input for a three-state multidimensional negative binomial HMM.
 * For each dimension `d`, there are two negative binomial distributions, low_d and high_d.
 * - ZERO state corresponds to zero emissions for all dimensions
 * - LOW state employs `low_d` emission for each dimension `d`
 * - HIGH state employs `high_d` emission for each dimension `d`
 *
 * @author Alexey Dievsky
 * @author Oleg Shpynov
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel> private constructor(
    fitInformation: Span1AnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int,
    saveExtendedInfo: Boolean = false
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, ZLH>(
    fitInformation, modelFitter, modelClass, ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L), fixedModelPath,
    threshold, maxIter, saveExtendedInfo
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

    companion object {

        const val SPAN_DEFAULT_BIN = 200
        const val SPAN_DEFAULT_FDR = 0.05
        const val SPAN_DEFAULT_GAP = 3

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
            multistartIter: Int = Fitter.MULTISTART_ITERATIONS,
            saveExtendedInfo: Boolean = false
        ): SpanPeakCallingExperiment<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.createFitInformation(
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
                    maxIter,
                    saveExtendedInfo
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
                    maxIter,
                    saveExtendedInfo
                )
            }
        }
    }
}

