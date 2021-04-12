package org.jetbrains.bio.experiments.fit.span234

import org.jetbrains.bio.experiments.fit.Span1AnalyzeFitInformation
import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * - ZERO state
 * - LOW state
 * - MEDIUM state
 * - HIGH state
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNBHMM4<Model : ClassificationModel> private constructor(
    fitInformation: Span1AnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, ZLMH>(
    fitInformation, modelFitter, modelClass, ZLMH.values(), NullHypothesis.of(ZLMH.Z, ZLMH.M, ZLMH.L), fixedModelPath,
    threshold, maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span-nbhmm4"

    companion object {

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
            multistartIter: Int = Fitter.MULTISTART_ITERATIONS
        ): SpanPeakCallingExperimentNBHMM4<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNBHMM4(
                fitInformation,
                when {
                    multistarts > 0 ->
                        NBHMM4.fitter().multiStarted(multistarts, multistartIter)
                    else ->
                        NBHMM4.fitter()
                },
                NBHMM4::class.java,
                fixedModelPath,
                threshold,
                maxIter
            )
        }
    }
}

enum class ZLMH {
    Z,  // ZERO
    L,  // LOW
    M,  // MEDIUM
    H;  // HIGH
}
