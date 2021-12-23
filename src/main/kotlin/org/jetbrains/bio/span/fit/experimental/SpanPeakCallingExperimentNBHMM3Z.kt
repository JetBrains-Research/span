package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.Span1AnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment
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
class SpanPeakCallingExperimentNBHMM3Z<Model : ClassificationModel> private constructor(
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

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.$MODEL_EXT"

    companion object {

        const val MODEL_EXT = "span-nbhmm3z"

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
        ): SpanPeakCallingExperimentNBHMM3Z<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNBHMM3Z(
                fitInformation,
                when {
                    multistarts > 0 ->
                        NBHMM3Z.fitter().multiStarted(multistarts, multistartIter)
                    else ->
                        NBHMM3Z.fitter()
                },
                NBHMM3Z::class.java,
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

class NBHMM3Z(nbMeans: DoubleArray, nbFailures: DoubleArray) : NBHMMZ(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NBHMM3Z> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM3Z = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM3Z {
                val (means, failures) = guess(preprocessed, 3, attempt)
                return NBHMM3Z(means, failures)
            }
        }
    }

}
