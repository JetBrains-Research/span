package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.statistics.hmm.FreeNBZHMM
import org.jetbrains.bio.statistics.Preprocessed
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
class SpanPeakCallingExperimentNB3ZHMM<Model : ClassificationModel> private constructor(
    fitInformation: SpanAnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int
) : SpanModelFitExperiment<Model, SpanAnalyzeFitInformation, ZLMH>(
    fitInformation, modelFitter, modelClass, ZLMH.values(), NullHypothesis.of(ZLMH.Z, ZLMH.M, ZLMH.L), fixedModelPath,
    threshold, maxIterations
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NB3Z_HMM.extension}"

    companion object {

        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            bin: Int,
            fragment: Fragment = AutoFragment,
            unique: Boolean = true,
            fixedModelPath: Path? = null,
            threshold: Double = Fitter.THRESHOLD,
            maxIterations: Int = Fitter.MAX_ITERATIONS,
            multistarts: Int = Fitter.MULTISTARTS,
            multistartIterations: Int = Fitter.MULTISTART_ITERATIONS
        ): SpanPeakCallingExperimentNB3ZHMM<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNB3ZHMM(
                fitInformation,
                when {
                    multistarts > 1 ->
                        NB3ZHMM.fitter().multiStarted(multistarts, multistartIterations)
                    else ->
                        NB3ZHMM.fitter()
                },
                NB3ZHMM::class.java,
                fixedModelPath,
                threshold,
                maxIterations
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

class NB3ZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : FreeNBZHMM(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NB3ZHMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int,
                attempt: Int
            ): NB3ZHMM = guess(listOf(preprocessed), title, threshold, maxIterations, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int,
                attempt: Int
            ): NB3ZHMM {
                val (means, failures) = guess(preprocessed, 3, attempt)
                return NB3ZHMM(means, failures)
            }
        }
    }

}
