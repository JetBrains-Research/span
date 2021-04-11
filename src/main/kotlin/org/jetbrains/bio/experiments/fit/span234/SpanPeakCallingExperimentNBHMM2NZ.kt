package org.jetbrains.bio.experiments.fit.span234

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.Span1AnalyzeFitInformation
import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * - LOW state
 * - HIGH state
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNBHMM2NZ<Model : ClassificationModel> private constructor(
    fitInformation: Span1AnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, LH>(
    fitInformation, modelFitter, modelClass, LH.values(), NullHypothesis.of(LH.L), fixedModelPath,
    threshold, maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

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
        ): SpanPeakCallingExperimentNBHMM2NZ<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported by the model" }
            return SpanPeakCallingExperimentNBHMM2NZ(
                fitInformation,
                when {
                    multistarts > 0 ->
                        NBHMM2NZ.fitter().multiStarted(multistarts, multistartIter)
                    else ->
                        NBHMM2NZ.fitter()
                },
                NBHMM2NZ::class.java,
                fixedModelPath,
                threshold,
                maxIter
            )
        }
    }
}

enum class LH {
    L,  // LOW
    H;  // HIGH
}


class NBHMM2NZ(nbMeans: DoubleArray, nbFailures: Double) : NBHMMNZ(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NBHMM2NZ> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM2NZ = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM2NZ {
                val (means, fs) = guess(preprocessed, 2, attempt)
                return NBHMM2NZ(means, fs)
            }
        }
    }

}
