package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.TrackAboutBooleanColumnType
import org.jetbrains.bio.genome.TrackAboutDoubleColumnType
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.genome.TrackAboutStringColumnType
import org.jetbrains.bio.span.fit.experimental.NB2HMM
import org.jetbrains.bio.span.fit.experimental.NB3HMM
import org.jetbrains.bio.span.fit.experimental.NB3ZHMM
import org.jetbrains.bio.span.fit.experimental.NB5ZHMM
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.span.statistics.mixture.NB2ZMixture
import org.jetbrains.bio.span.statistics.mixture.NegBin2ZeroRegressionMixture
import org.jetbrains.bio.span.statistics.mixture.PoissonRegression2Mixture
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * Contains the results of a Span-like model-fitting experiment.
 *
 * @property fitInfo The [SpanFitInformation] instance that describes the experiment input.
 * @property model The [ClassificationModel] that was fitted during the experiment.
 * @property logNullMemberships The chromosome-wise dataframes of log null probabilities, i.e.
 * the log probability of each observation under the null hypothesis. Each dataframe should at least contain
 * a column of floats or doubles labelled [SpanModelFitExperiment.NULL].
 * @property statesDataFrameMap extended information which contains additional states per position mapping
 */
open class SpanFitResults(
    val fitInfo: SpanFitInformation,
    val model: ClassificationModel,
    val logNullMemberships: Map<String, DataFrame>,
    val statesDataFrameMap: Map<String, DataFrame>?
) {
    companion object {
        internal val LOG = LoggerFactory.getLogger(SpanFitResults::class.java)

        val CT_MODEL_FILE = TrackAboutStringColumnType("Model")
        val CT_MODEL_TYPE = TrackAboutStringColumnType("Model type")
        val CT_SIGNAL_MEAN = TrackAboutDoubleColumnType("Signal mean")
        val CT_NOISE_MEAN = TrackAboutDoubleColumnType("Noise mean")
        val CT_SIGNAL_TO_NOISE = TrackAboutDoubleColumnType("Signal to noise")
        val CT_OUT_OF_SNR_UP = TrackAboutBooleanColumnType("Out of signal-to-noise range up")
        val CT_OUT_OF_SNR_DOWN = TrackAboutBooleanColumnType("Out of signal-to-noise range down")
    }

    /**
     * @return Information about fit results including model and other parameters
     */
    open fun modelInformation(modelPath: Path): List<TrackAboutMetricValue<*>> {
        return when (model) {
            is NB2ZHMM -> {
                val outOfSnrHitUp = model.outOfSignalToNoiseRatioRangeUp
                val outOfSnrHitDown = model.outOfSignalToNoiseRatioRangeDown
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                listOf(
                    CT_MODEL_FILE to modelPath,
                    CT_MODEL_TYPE to SpanModelType.NB2Z_HMM.description,
                    CT_SIGNAL_MEAN to signalMean,
                    CT_NOISE_MEAN to noiseMean,
                    CT_SIGNAL_TO_NOISE to ((signalMean + 1e-10) / (noiseMean + 1e-10)),
                    CT_OUT_OF_SNR_UP to outOfSnrHitUp,
                    CT_OUT_OF_SNR_DOWN to outOfSnrHitDown,
                )
            }

            is NB2ZMixture -> {
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                listOf(
                    CT_MODEL_TYPE to SpanModelType.NB2Z_MIXTURE.description,
                    CT_SIGNAL_MEAN to signalMean,
                    CT_NOISE_MEAN to noiseMean,
                    CT_SIGNAL_TO_NOISE to ((signalMean + 1e-10) / (noiseMean + 1e-10)),
                )
            }

            is PoissonRegression2Mixture -> listOf(
                CT_MODEL_TYPE to SpanModelType.POISSON_REGRESSION_MIXTURE.description,
                CT_SIGNAL_TO_NOISE to model.signalToNoise,
            )

            is NegBin2ZeroRegressionMixture -> listOf(
                CT_MODEL_TYPE to SpanModelType.NEGBIN_REGRESSION_MIXTURE.description,
                CT_SIGNAL_TO_NOISE to model.signalToNoise,
            )

            is NB3ZHMM -> listOf(CT_MODEL_TYPE to SpanModelType.NB3Z_HMM.description)
            is NB5ZHMM -> listOf(CT_MODEL_TYPE to SpanModelType.NB5Z_HMM.description)
            is NB2HMM -> listOf(CT_MODEL_TYPE to SpanModelType.NB2_HMM.description)
            is NB3HMM -> listOf(CT_MODEL_TYPE to SpanModelType.NB3_HMM.description)
            else -> emptyList()
        }
    }
}
