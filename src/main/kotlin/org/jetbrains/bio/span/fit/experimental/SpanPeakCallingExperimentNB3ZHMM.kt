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
import org.jetbrains.bio.span.statistics.hmm.FreeNBZHMM
import org.jetbrains.bio.span.statistics.util.NegBinUtil.guessByData
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
    maxIterations: Int,
    saveExtendedInfo: Boolean,
    keepModelFile: Boolean
) : SpanModelFitExperiment<Model, SpanAnalyzeFitInformation, ZLMH>(
    fitInformation, modelFitter, modelClass,
    ZLMH.values(), NullHypothesis.of(ZLMH.Z, ZLMH.M, ZLMH.L),
    threshold, maxIterations,
    fixedModelPath, saveExtendedInfo, keepModelFile
) {

    override val defaultModelPath: Path =
        experimentPath / "${fitInformation.id}.${SpanModelType.NB3Z_HMM.extension}"

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
        ): SpanPeakCallingExperimentNB3ZHMM<out ClassificationModel> {
            require(paths.isNotEmpty()) { "No data" }
            val fitInformation = SpanAnalyzeFitInformation.createFitInformation(
                genomeQuery, paths, explicitFormat, MultiLabels.generate(SPAN_TRACK_PREFIX, paths.size).toList(),
                fragment, unique, bin
            )
            require(paths.size == 1) { "Multiple replicates are not supported" }
            return SpanPeakCallingExperimentNB3ZHMM(
                fitInformation,
                NB3ZHMM.fitter(),
                NB3ZHMM::class.java,
                fixedModelPath,
                threshold,
                maxIterations,
                saveExtendedInfo,
                keepModelFile
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
                maxIterations: Int
            ): NB3ZHMM = guess(listOf(preprocessed), title, threshold, maxIterations)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB3ZHMM {
                val guess = guessByData(positiveCoverage(preprocessed), 3)
                return NB3ZHMM(guess.means, guess.failures)
            }
        }
    }

}
