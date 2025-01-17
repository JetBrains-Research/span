package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.span.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Span `analyze --type nbhmm` invocation.
 *
 * For each treatment-control pair, we compute binned normalized coverage.
 * These coverages are used as the input for a three-state multidimensional negative binomial HMM.
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
    fitInformation: SpanAnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepCacheFiles: Boolean
) : SpanModelFitExperiment<Model, SpanAnalyzeFitInformation, ZLH>(
    fitInformation, modelFitter, modelClass, ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepCacheFiles
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NB2Z_HMM.extension}"

    companion object {

        /**
         * Technical prefix used for generating track labels.
         */
        const val SPAN_TRACK_PREFIX = "track"

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
            hmmEstimateSNR: Double = SPAN_DEFAULT_HMM_ESTIMATE_SNR,
            hmmLow: Double = SPAN_DEFAULT_HMM_LOW_THRESHOLD,
            fixedModelPath: Path? = null,
            threshold: Double = SPAN_DEFAULT_FIT_THRESHOLD,
            maxIterations: Int = SPAN_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo: Boolean = false,
            keepCacheFiles: Boolean = false
        ): SpanPeakCallingExperiment<out ClassificationModel> {
            require(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, MultiLabels.generate(SPAN_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            return if (paths.size == 1) {
                SpanPeakCallingExperiment(
                    fitInformation,
                    NB2ZHMM.fitter(hmmEstimateSNR, hmmLow),
                    NB2ZHMM::class.java,
                    fixedModelPath,
                    threshold, maxIterations,
                    saveExtendedInfo,
                    keepCacheFiles
                )
            } else {
                SpanPeakCallingExperiment(
                    fitInformation,
                    ConstrainedNBZHMM.fitter(paths.size),
                    ConstrainedNBZHMM::class.java,
                    fixedModelPath,
                    threshold, maxIterations,
                    saveExtendedInfo,
                    keepCacheFiles
                )
            }
        }
    }
}

