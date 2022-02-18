package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation.Companion.createFitInformation
import org.jetbrains.bio.span.statistics.mixture.NB2ZMixture
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path


/**
 * - LOW state
 * - MEDIUM state
 * - HIGH state
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNB2ZMixture private constructor(
    fitInformation: SpanAnalyzeFitInformation,
    fitter: Fitter<NB2ZMixture>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean = true
) : SpanModelFitExperiment<NB2ZMixture, SpanAnalyzeFitInformation, ZLH>(
    fitInformation,
    fitter,
    NB2ZMixture::class.java,
    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
    fixedModelPath, threshold, maxIterations, saveExtendedInfo
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NB2Z_MIXTURE.extension}"

    companion object {

        /**
         * Contains a check that a single treatment-control pair was provided.
         */
        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            fragment: Fragment,
            binSize: Int,
            unique: Boolean,
            fixedModelPath: Path?,
            threshold: Double = Fitter.THRESHOLD,
            maxIterations: Int = Fitter.MAX_ITERATIONS,
            multistarts: Int = Fitter.MULTISTARTS,
            multistartIterations: Int = Fitter.MULTISTART_ITERATIONS,
            saveExtendedInfo: Boolean = false
        ): SpanPeakCallingExperimentNB2ZMixture {
            require(paths.size == 1) { "Mixture currently accepts a single data track." }
            val fitInformation = createFitInformation(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, binSize
            )
            return SpanPeakCallingExperimentNB2ZMixture(
                fitInformation,
                when {
                    multistarts > 1 ->
                        NB2ZMixture.fitter().multiStarted(multistarts, multistartIterations)
                    else ->
                        NB2ZMixture.fitter()
                },
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo
            )
        }
    }
}
