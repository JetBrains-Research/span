package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_LOW_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_PRIORS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_TRANSITIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_NB_VAR_MEAN_MULTIPLIER
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.span.statistics.util.NegBinUtil
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution.Companion.estimateFailuresUsingMoments
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import kotlin.math.max

/**
 * A zero-inflated HMM with univariate Negative Binomial emissions and constraints.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 25/05/15
 */
class NB2ZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) :
    FreeNBZHMM(nbMeans, nbFailures,
        priorProbabilities = SPAN_HMM_PRIORS,
        transitionProbabilities = F64Array.invoke(3, 3) { i, j -> SPAN_HMM_TRANSITIONS[i][j] }
    ) {

    // Indicator: signal-to-noise ratio is smaller than estimation
    var outOfSignalToNoiseRatioRangeDown: Boolean = false
    // Indicator: low state mean value is smaller than threshold
    var outOfLowerNoise: Boolean = false

    /**
     * Keep model signal-to-noise ratio in the normal range
     */
    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        super.updateParameters(df, gammas)

        // You wanna keep em separated!
        for (d in 0 until numDimensions) {
            val lowState = getEmissionScheme(1, d) as NegBinEmissionScheme
            val highState = getEmissionScheme(2, d) as NegBinEmissionScheme

            // Need to update transients in case of any change
            var updated = false

            // This check is required to prevent low state go too close to 0, causing too broad peaks
            val lowMeanThreshold = guess.lowMean * SPAN_HMM_LOW_THRESHOLD
            if (lowState.mean < lowMeanThreshold) {
                LOG.warn("Low state mean ${lowState.mean} < ${lowMeanThreshold}, fixing...")
                outOfLowerNoise = true
                lowState.mean = lowMeanThreshold
                lowState.failures = estimateFailuresUsingMoments(
                    lowState.mean,
                    max(lowState.mean * SPAN_HMM_NB_VAR_MEAN_MULTIPLIER, lowState.variance)
                )
                updated = true
            }

            // This check is required mostly for narrow marks with low SNR
            val snr = (highState.mean + 1e-10) / (lowState.mean + 1e-10)
            if (snr < guess.signalToNoise) {
                LOG.warn("Signal-to-noise ratio $snr < ${guess.signalToNoise}, fixing...")
                outOfSignalToNoiseRatioRangeDown = true
                highState.mean = lowState.mean * guess.signalToNoise
                highState.failures = estimateFailuresUsingMoments(
                    highState.mean,
                    max(highState.mean * SPAN_HMM_NB_VAR_MEAN_MULTIPLIER, highState.variance)
                )
                updated = true
            }

            if (updated) {
                highState.updateTransients()
                lowState.updateTransients()
            }
        }
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 3

        // This will be updated during guess step
        lateinit var guess: NegBinUtil.Guess

        private val LOG: Logger = LoggerFactory.getLogger(NB2ZHMM::class.java)

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
                guess = guess(preprocessed, 2)
                return NB2ZHMM(guess.means, guess.failures)
            }
        }
    }
}
