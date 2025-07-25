package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation.Companion.createFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.SPAN_TRACK_PREFIX
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
    saveExtendedInfo: Boolean,
    keepModelFile: Boolean
) : SpanModelFitExperiment<NB2ZMixture, SpanAnalyzeFitInformation, ZLH>(
    fitInformation, fitter, NB2ZMixture::class.java,
    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepModelFile
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
            explicitFormat: ReadsFormat?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int,
            fixedModelPath: Path?,
            threshold: Double = SPAN_DEFAULT_FIT_THRESHOLD,
            maxIterations: Int = SPAN_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo: Boolean = false,
            keepModelFile: Boolean = false
        ): SpanPeakCallingExperimentNB2ZMixture {
            require(paths.size == 1) { "Mixture currently accepts a single data track." }
            val fitInformation = createFitInformation(
                genomeQuery, paths, explicitFormat, MultiLabels.generate(SPAN_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, binSize
            )
            return SpanPeakCallingExperimentNB2ZMixture(
                fitInformation,
                NB2ZMixture.fitter(),
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo,
                keepModelFile
            )
        }
    }
}
