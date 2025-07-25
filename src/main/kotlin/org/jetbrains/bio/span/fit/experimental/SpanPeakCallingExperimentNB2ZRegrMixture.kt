package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.fit.SpanModelType
import org.jetbrains.bio.span.fit.ZLH
import org.jetbrains.bio.span.statistics.mixture.NegBin2ZeroRegressionMixture
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Span `analyze-experimental --type nbrm` invocation.
 *
 * Currently, supports only a single treatment track.
 *
 * We compute binned coverage for the treatment track and use it as the response vector.
 *
 * We also compute the covariates:
 * - "GC" and "GC2" are binned mean GC content and its square
 * - "input" is the binned control track coverage, if supplied
 * - "mapability" is the binned mean mapability, if supplied
 *
 * These data are used as the input for a three-state negative binomial regression mixture.
 * - ZERO state corresponds to zero emission
 * - LOW state employs a negative binomial GLM with the covariates listed above
 * - HIGH state employs another negative binomial GLM with the covariates listed above
 */
class SpanPeakCallingExperimentNB2ZRegrMixture private constructor(
    fitInformation: SpanRegrMixtureAnalyzeFitInformation,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepModelFile: Boolean
) : SpanModelFitExperiment<NegBin2ZeroRegressionMixture, SpanRegrMixtureAnalyzeFitInformation, ZLH>(
    fitInformation, NegBin2ZeroRegressionMixture.fitter(), NegBin2ZeroRegressionMixture::class.java,
    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepModelFile
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NEGBIN_REGRESSION_MIXTURE.extension}"

    companion object {

        /**
         * Contains a check that a single treatment-control pair was provided.
         */
        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            explicitFormat: ReadsFormat?,
            mapabilityPath: Path?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int,
            fixedModelPath: Path?,
            threshold: Double,
            maxIterations: Int,
            saveExtendedInfo: Boolean,
            keepModelFile: Boolean
        ): SpanPeakCallingExperimentNB2ZRegrMixture {
            require(paths.size == 1) { "Negative binomial regression mixture currently accepts a single data track." }
            val fitInformation = SpanRegrMixtureAnalyzeFitInformation(
                genomeQuery, paths.single(), explicitFormat, mapabilityPath, fragment, unique, binSize
            )
            return SpanPeakCallingExperimentNB2ZRegrMixture(
                fitInformation,
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo,
                keepModelFile
            )
        }
    }
}

