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
 * - MEDIUM 1 state
 * - MEDIUM 2 state
 * - HIGH state
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNBHMM4Z<Model : ClassificationModel> private constructor(
    fitInformation: Span1AnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIter: Int
) : SpanModelFitExperiment<Model, Span1AnalyzeFitInformation, ZLM2H>(
    fitInformation,
    modelFitter,
    modelClass,
    ZLM2H.values(),
    NullHypothesis.of(ZLM2H.Z, ZLM2H.M1, ZLM2H.M2, ZLM2H.L),
    fixedModelPath,
    threshold,
    maxIter
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.$MODEL_EXT"

    companion object {
        const val MODEL_EXT = "span-nbhmm4z"

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
        ): SpanPeakCallingExperimentNBHMM4Z<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = Span1AnalyzeFitInformation.effective(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNBHMM4Z(
                fitInformation,
                when {
                    multistarts > 0 ->
                        NBHMM4Z.fitter().multiStarted(multistarts, multistartIter)
                    else ->
                        NBHMM4Z.fitter()
                },
                NBHMM4Z::class.java,
                fixedModelPath,
                threshold,
                maxIter
            )
        }
    }
}

enum class ZLM2H {
    Z,  // ZERO
    L,  // LOW
    M1,  // MEDIUM1
    M2,  // MEDIUM2
    H;  // HIGH
}

class NBHMM4Z(nbMeans: DoubleArray, nbFailures: DoubleArray) : NBHMMZ(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NBHMM4Z> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM4Z = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM4Z {
                val (means, failures) = guess(preprocessed, 4, attempt)
                return NBHMM4Z(means, failures)
            }
        }
    }

}
