package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_THRESHOLD
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
 * - (LOW) S1 - (HIGH) S5 states
 *
 * @author Oleg Shpynov
 */
class SpanPeakCallingExperimentNB5ZHMM<Model : ClassificationModel> private constructor(
    fitInformation: SpanAnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepModelFile: Boolean
) : SpanModelFitExperiment<Model, SpanAnalyzeFitInformation, Z5>(
    fitInformation, modelFitter, modelClass,
    Z5.values(), NullHypothesis.of(Z5.Z, Z5.S1, Z5.S2, Z5.S3, Z5.S4),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepModelFile
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
            threshold: Double = SPAN_DEFAULT_FIT_THRESHOLD,
            maxIterations: Int = SPAN_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo: Boolean = false,
            keepModelFile: Boolean = false
        ): SpanPeakCallingExperimentNB5ZHMM<out ClassificationModel> {
            require(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.SPAN_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNB5ZHMM(
                fitInformation,
                NB5ZHMM.fitter(),
                NB5ZHMM::class.java,
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo,
                keepModelFile
            )
        }
    }
}

enum class Z5 {
    Z,  // ZERO
    S1,  // LOW
    S2,
    S3,
    S4,
    S5;  // HIGH
}

class NB5ZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : FreeNBZHMM(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NB5ZHMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB5ZHMM = guess(listOf(preprocessed), title, threshold, maxIterations)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB5ZHMM {
                val guess = guess(preprocessed, 5)
                return NB5ZHMM(guess.means, guess.failures)
            }
        }
    }

}
