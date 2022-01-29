package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation.Companion.createFitInformation
import org.jetbrains.bio.span.statistics.mixture.NegBinMixture
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
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
class SpanPeakCallingExperimentNB2ZM private constructor(
    fitInformation: SpanAnalyzeFitInformation,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int,
    saveExtendedInfo: Boolean = true
) : SpanModelFitExperiment<NegBinMixture, SpanAnalyzeFitInformation, ZLH>(
    fitInformation,
    NegBinMixture.fitter(), NegBinMixture::class.java,
    ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
    fixedModelPath, threshold, maxIter, saveExtendedInfo
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
            threshold: Double,
            maxIter: Int,
            saveExtendedInfo: Boolean = false
        ): SpanPeakCallingExperimentNB2ZM {
            check(paths.size == 1) { "Mixture currently accepts a single data track." }
            val fitInformation = createFitInformation(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, binSize
            )
            return SpanPeakCallingExperimentNB2ZM(fitInformation, fixedModelPath, threshold, maxIter, saveExtendedInfo)
        }
    }
}
