package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.*
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
class SpanPeakCallingExperimentNB2HMM<Model : ClassificationModel> private constructor(
    fitInformation: SpanAnalyzeFitInformation,
    modelFitter: Fitter<Model>,
    modelClass: Class<Model>,
    fixedModelPath: Path?,
    threshold: Double,
    maxIterations: Int,
    saveExtendedInfo: Boolean
) : SpanModelFitExperiment<Model, SpanAnalyzeFitInformation, LH>(
    fitInformation, modelFitter, modelClass, LH.values(), NullHypothesis.of(LH.L), fixedModelPath,
    threshold, maxIterations, saveExtendedInfo = saveExtendedInfo
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NB2_HMM.extension}"

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
            multistartIterations: Int = Fitter.MULTISTART_ITERATIONS,
            saveExtendedInfo: Boolean
        ): SpanPeakCallingExperimentNB2HMM<out ClassificationModel> {
            check(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, MultiLabels.generate(SpanPeakCallingExperiment.TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported by the model" }
            return SpanPeakCallingExperimentNB2HMM(
                fitInformation,
                when {
                    multistarts > 1 ->
                        NB2HMM.fitter().multiStarted(multistarts, multistartIterations)
                    else ->
                        NB2HMM.fitter()
                },
                NB2HMM::class.java,
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo
            )
        }
    }
}

enum class LH {
    L,  // LOW
    H;  // HIGH
}


class NB2HMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : FreeNBHMM(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1


        fun fitter() = object : Fitter<NB2HMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int,
                attempt: Int
            ): NB2HMM = guess(listOf(preprocessed), title, threshold, maxIterations, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int,
                attempt: Int
            ): NB2HMM {
                val (means, fs) = guess(preprocessed, 2, attempt)
                return NB2HMM(means, fs)
            }
        }
    }

}
