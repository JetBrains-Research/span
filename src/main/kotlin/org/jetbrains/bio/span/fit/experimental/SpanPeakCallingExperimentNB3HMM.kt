package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.SPAN_TRACK_PREFIX
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
class SpanPeakCallingExperimentNB3HMM<Model : ClassificationModel> private constructor(
    fitInformation: SpanAnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepModelFile: Boolean
) : SpanModelFitExperiment<Model, SpanAnalyzeFitInformation, LMH>(
    fitInformation, modelFitter, modelClass,
    LMH.values(), NullHypothesis.of(LMH.L, LMH.M),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepModelFile
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NB3_HMM.extension}"

    companion object {

        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            explicitFormat: ReadsFormat?,
            fragment: Fragment = AutoFragment,
            unique: Boolean = true,
            bin: Int,
            fixedModelPath: Path? = null,
            threshold: Double = SPAN_DEFAULT_FIT_THRESHOLD,
            maxIterations: Int = SPAN_DEFAULT_FIT_MAX_ITERATIONS,
            saveExtendedInfo: Boolean = false,
            keepModelFile: Boolean = false
        ): SpanPeakCallingExperimentNB3HMM<out ClassificationModel> {
            require(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, explicitFormat, MultiLabels.generate(SPAN_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNB3HMM(
                fitInformation,
                NB3HMM.fitter(),
                NB3HMM::class.java,
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo,
                keepModelFile
            )
        }
    }
}

enum class LMH {
    L,  // LOW
    M,  // MEDIUM
    H;  // HIGH
}

class NB3HMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : FreeNBHMM(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NB3HMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB3HMM = guess(listOf(preprocessed), title, threshold, maxIterations)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB3HMM {
                val (means, fs) = guess(preprocessed, 3)
                return NB3HMM(means, fs)
            }
        }
    }

}
