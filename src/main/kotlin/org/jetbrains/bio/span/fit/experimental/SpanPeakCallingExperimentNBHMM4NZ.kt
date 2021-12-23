package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.Span1AnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.TRACK_PREFIX
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * - ALMOST_ZERO state
 * - LOW state
 * - MEDIUM state
 * - HIGH state
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNBHMM4NZ<Model : ClassificationModel> private constructor(
    fitInformation: Span1AnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, AZLMH>(
    fitInformation,
    modelFitter,
    modelClass,
    AZLMH.values(),
    NullHypothesis.of(AZLMH.A, AZLMH.L, AZLMH.M),
    fixedModelPath,
    threshold,
    maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.$MODEL_EXT"

    companion object {

        const val MODEL_EXT = "span-nbhmm4nz"

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
        ): SpanPeakCallingExperimentNBHMM4NZ<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNBHMM4NZ(
                fitInformation,
                when {
                    multistarts > 0 ->
                        NBHMM4NZ.fitter().multiStarted(multistarts, multistartIter)
                    else ->
                        NBHMM4NZ.fitter()
                },
                NBHMM4NZ::class.java,
                fixedModelPath,
                threshold,
                maxIter
            )
        }
    }
}

enum class AZLMH {
    A,  // ALMOST ZERO
    L,  // LOW
    M,  // MEDIUM
    H;  // HIGH
}

class NBHMM4NZ(nbMeans: DoubleArray, nbFailures: Double) : NBHMMNZ(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NBHMM4NZ> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM4NZ = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM4NZ {
                val (means, fs) = guess(preprocessed, 4, attempt)
                return NBHMM4NZ(means, fs)
            }
        }
    }

}
