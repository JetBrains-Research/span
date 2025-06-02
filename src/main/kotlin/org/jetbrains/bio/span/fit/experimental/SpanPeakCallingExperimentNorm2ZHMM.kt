package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.span.statistics.hmm.Norm2ZHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

class SpanPeakCallingExperimentNorm2ZHMM<Model : ClassificationModel> private constructor(
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
        experimentPath / "${fitInformation.id}.${SpanModelType.NORM2Z_HMM.extension}"

    companion object {

        const val SPAN_TRACK_PREFIX = "track"

        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            explicitFormat: ReadsFormat?,
            fragment: Fragment = AutoFragment,
            unique: Boolean = true,
            bin: Int,
            hmmEstimateSNR: Double = SPAN_DEFAULT_HMM_ESTIMATE_SNR,
            hmmLow: Double = SPAN_DEFAULT_HMM_LOW_THRESHOLD,
            fixedModelPath: Path? = null,
            threshold: Double = SPAN_DEFAULT_FIT_THRESHOLD,
            maxIterations: Int = SPAN_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo: Boolean = false,
            keepCacheFiles: Boolean = false
        ): SpanPeakCallingExperimentNorm2ZHMM<out ClassificationModel> {
            require(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, explicitFormat, MultiLabels.generate(SPAN_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            return if (paths.size == 1) {
                SpanPeakCallingExperimentNorm2ZHMM(
                    fitInformation,
                    Norm2ZHMM.fitter(hmmEstimateSNR, hmmLow),
                    Norm2ZHMM::class.java,
                    fixedModelPath,
                    threshold, maxIterations,
                    saveExtendedInfo,
                    keepCacheFiles
                )
            } else {
                throw UnsupportedOperationException("Multiple replicates are not supported by the model")
            }
        }
    }
}

