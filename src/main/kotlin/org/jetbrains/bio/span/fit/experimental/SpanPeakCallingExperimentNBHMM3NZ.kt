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
 * - LOW state
 * - MEDIUM state
 * - HIGH state
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNBHMM3NZ<Model : ClassificationModel> private constructor(
    fitInformation: Span1AnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, LMH>(
    fitInformation, modelFitter, modelClass, LMH.values(), NullHypothesis.of(LMH.L, LMH.M), fixedModelPath,
    threshold, maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.$MODEL_EXT"

    companion object {

        const val MODEL_EXT = "span-nbhmm3nz"

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
        ): SpanPeakCallingExperimentNBHMM3NZ<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNBHMM3NZ(
                fitInformation,
                when {
                    multistarts > 0 ->
                        NBHMM3NZ.fitter().multiStarted(multistarts, multistartIter)
                    else ->
                        NBHMM3NZ.fitter()
                },
                NBHMM3NZ::class.java,
                fixedModelPath,
                threshold,
                maxIter
            )
        }
    }
}

enum class LMH {
    L,  // LOW
    M,  // MEDIUM
    H;  // HIGH
}

class NBHMM3NZ(nbMeans: DoubleArray, nbFailures: Double) : NBHMMNZ(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NBHMM3NZ> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM3NZ = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM3NZ {
                val (means, fs) = guess(preprocessed, 3, attempt)
                return NBHMM3NZ(means, fs)
            }
        }
    }

}
