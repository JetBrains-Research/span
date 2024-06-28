package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_NB2ZHMM_PRIORS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_NB2ZHMM_TRANSITIONS_
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MAX_SIGNAL_TO_NOISE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MIN_SIGNAL_TO_NOISE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SIGNAL_TO_NOISE_PUSHBACK
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import kotlin.math.sqrt

/**
 * A zero-inflated HMM with univariate Negative Binomial emissions.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 25/05/15
 */
class NB2ZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) :
    FreeNBZHMM(nbMeans, nbFailures,
        priorProbabilities = SPAN_DEFAULT_NB2ZHMM_PRIORS,
        transitionProbabilities = F64Array.invoke(3, 3) { i, j -> SPAN_DEFAULT_NB2ZHMM_TRANSITIONS_[i][j] }
    ) {

    // Indicator value to save that model fit was corrected during fitting
    var outOfSignalToNoiseRatioRangeUp: Boolean = false
    var outOfSignalToNoiseRatioRangeDown: Boolean = false

    private var prevDistributions: List<Pair<NegBinEmissionScheme, NegBinEmissionScheme>> = emptyList()

    /**
     * Keep model signal-to-noise ratio in the normal range
     */
    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        super.updateParameters(df, gammas)

        // You wanna keep em separated!
        for (d in 0 until numDimensions) {
            val lowState = getEmissionScheme(1, d) as NegBinEmissionScheme
            val highState = getEmissionScheme(2, d) as NegBinEmissionScheme
            val snr = highState.mean / lowState.mean
            LOG.debug("Signal-to-noise ratio $snr")
            if (snr !in SPAN_MIN_SIGNAL_TO_NOISE..SPAN_MAX_SIGNAL_TO_NOISE) {
                LOG.warn(
                    "Signal-to-noise ratio not in $SPAN_MIN_SIGNAL_TO_NOISE-$SPAN_MAX_SIGNAL_TO_NOISE, fixing..."
                )

                val gMean = if (prevDistributions.isNotEmpty())
                    sqrt(prevDistributions[d].first.mean * prevDistributions[d].second.mean)
                else
                    sqrt(lowState.mean * highState.mean)
                val snrRange = SPAN_MAX_SIGNAL_TO_NOISE - SPAN_MIN_SIGNAL_TO_NOISE
                val targetSnr: Double
                if (snr < SPAN_MIN_SIGNAL_TO_NOISE) {
                    outOfSignalToNoiseRatioRangeDown = true
                    targetSnr = SPAN_MIN_SIGNAL_TO_NOISE + SPAN_SIGNAL_TO_NOISE_PUSHBACK * snrRange
                } else {
                    outOfSignalToNoiseRatioRangeUp = true
                    targetSnr = SPAN_MAX_SIGNAL_TO_NOISE - SPAN_SIGNAL_TO_NOISE_PUSHBACK * snrRange
                }

                // We update mean values
                lowState.mean = gMean / sqrt(targetSnr)
                highState.mean = gMean * sqrt(targetSnr)

                // But keep failures shape from the previous step
                if (prevDistributions.isNotEmpty()) {
                    lowState.failures = prevDistributions[d].first.failures
                    highState.failures = prevDistributions[d].second.failures
                }

                highState.updateTransients()
                lowState.updateTransients()
            }
        }
        // Save current distributions
        prevDistributions = (0 until numDimensions).map { d ->
            val lowState = getEmissionScheme(1, d) as NegBinEmissionScheme
            val highState = getEmissionScheme(2, d) as NegBinEmissionScheme
            lowState to highState
        }.toList()
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 3

        private val LOG: Logger = org.slf4j.LoggerFactory.getLogger(NB2ZHMM::class.java)

        fun fitter() = object : Fitter<NB2ZHMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB2ZHMM = guess(listOf(preprocessed), title, threshold, maxIterations)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB2ZHMM {
                val (means, failures) = guess(preprocessed, 2)
                return NB2ZHMM(means, failures)
            }
        }
    }
}
