package org.jetbrains.bio.experiments.fit.experimental

import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.mixture.PoissonRegressionMixture
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Span `analyze --type prm` invocation.
 *
 * Currently supports only a single treatment track.
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
class Span2PeakCallingExperiment private constructor(
        fitInformation: Span2FitInformation,
        fixedModelPath: Path?,
        threshold: Double,
        maxIter: Int
) : SpanModelFitExperiment<PoissonRegressionMixture, Span2FitInformation, ZLH>(
        fitInformation,
        PoissonRegressionMixture.fitter(), PoissonRegressionMixture::class.java,
        ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
        fixedModelPath,
        threshold, maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span2"

    companion object {

        /**
         * Contains a check that a single treatment-control pair was provided.
         */
        fun getExperiment(
                genomeQuery: GenomeQuery,
                data: List<SpanDataPaths>,
                mapabilityPath: Path?,
                fragment: Fragment,
                binSize: Int,
                unique: Boolean,
                fixedModelPath: Path?,
                threshold: Double,
                maxIter: Int
        ): Span2PeakCallingExperiment {
            check(data.size == 1) { "Poisson regression mixture currently accepts a single data track." }
            val fitInformation = Span2FitInformation(
                    genomeQuery, data.single(), mapabilityPath, fragment, unique, binSize
            )
            return Span2PeakCallingExperiment(fitInformation, fixedModelPath, threshold, maxIter)
        }
    }
}
