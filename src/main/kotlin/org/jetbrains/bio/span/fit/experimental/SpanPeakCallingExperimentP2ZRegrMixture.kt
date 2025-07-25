package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.fit.SpanModelType
import org.jetbrains.bio.span.fit.ZLH
import org.jetbrains.bio.span.statistics.mixture.PoissonRegression2Mixture
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Span `analyze --type prm` invocation.
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
 * These data are used as the input for a three-state Poisson regression mixture.
 * - ZERO state corresponds to zero emission
 * - LOW state employs a Poisson GLM with the covariates listed above
 * - HIGH state employs another Poisson GLM with the covariates listed above
 */
class SpanPeakCallingExperimentP2ZRegrMixture private constructor(
    fitInformation: SpanRegrMixtureAnalyzeFitInformation,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepModelFile: Boolean
) : SpanModelFitExperiment<PoissonRegression2Mixture, SpanRegrMixtureAnalyzeFitInformation, ZLH>(
    fitInformation, PoissonRegression2Mixture.fitter(), PoissonRegression2Mixture::class.java,
    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepModelFile
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.POISSON_REGRESSION_MIXTURE.extension}"

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
            keepModelFile: Boolean = false
        ): SpanPeakCallingExperimentP2ZRegrMixture {
            require(paths.size == 1) { "Poisson regression mixture currently accepts a single data track." }
            val fitInformation = SpanRegrMixtureAnalyzeFitInformation(
                genomeQuery, paths.single(), explicitFormat, mapabilityPath, fragment, unique, binSize
            )
            return SpanPeakCallingExperimentP2ZRegrMixture(
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
