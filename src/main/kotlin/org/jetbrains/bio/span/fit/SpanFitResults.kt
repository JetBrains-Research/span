package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.TrackAboutDoubleColumnType
import org.jetbrains.bio.genome.TrackAboutMetricValue
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.span.statistics.mixture.NegBinMixture
import org.jetbrains.bio.span.statistics.mixture.NegBinRegressionMixture
import org.jetbrains.bio.span.statistics.mixture.PoissonRegressionMixture
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.slf4j.LoggerFactory

/**
 * Contains the results of a Span-like model-fitting experiment.
 *
 * @property fitInfo The [SpanFitInformation] instance that describes the experiment input.
 * @property model The [ClassificationModel] that was fitted during the experiment.
 * @property logNullMemberships The chromosome-wise dataframes of log null probabilities, i.e.
 * the log probability of each observation under the null hypothesis. Each dataframe should at least contain
 * a column of floats or doubles labelled [SpanModelFitExperiment.NULL].
 */
open class SpanFitResults(
    val fitInfo: SpanFitInformation,
    val model: ClassificationModel,
    val logNullMemberships: Map<String, DataFrame>
) {
    companion object {
        internal val LOG = LoggerFactory.getLogger(SpanFitResults::class.java)

        val CT_SIGNAL_MEAN = TrackAboutDoubleColumnType("Signal mean")
        val CT_NOISE_MEAN = TrackAboutDoubleColumnType("Noise mean")
        val CT_SIGNAL_TO_NOISE = TrackAboutDoubleColumnType("Signal to noise")
    }

    /**
     * @return Information about fit results including model and other parameters
     */
    open fun about(): List<TrackAboutMetricValue<*>> {
        return when (model) {
            is NB2ZHMM -> {
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                listOf(
                    CT_SIGNAL_MEAN to signalMean,
                    CT_NOISE_MEAN to noiseMean,
                    CT_SIGNAL_TO_NOISE to ((signalMean + 1e-10) / (noiseMean + 1e-10))
                )
            }
            is NegBinMixture -> {
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                listOf(
                    CT_SIGNAL_MEAN to signalMean,
                    CT_NOISE_MEAN to noiseMean,
                    CT_SIGNAL_TO_NOISE to ((signalMean + 1e-10) / (noiseMean + 1e-10))
                )
            }
            is PoissonRegressionMixture -> listOf(CT_SIGNAL_TO_NOISE to model.signalToNoise)
            is NegBinRegressionMixture -> listOf(CT_SIGNAL_TO_NOISE to model.signalToNoise)
            else -> emptyList()
        }
    }
}
