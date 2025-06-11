package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_NB_VAR_MEAN_MULTIPLIER
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_PRIORS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_TRANSITIONS
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.span.statistics.util.NegBinUtil
import org.jetbrains.bio.span.statistics.util.NegBinUtil.guessByData
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution.Companion.estimateFailuresUsingMoments
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import kotlin.math.max
import kotlin.math.sqrt

/**
 * A zero-inflated constrained HMM with univariate Negative Binomial emissions and constraints.
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
            var updatedLow = false
            var updatedHigh = false

            val snrPrevious = highState.mean / lowState.mean
            val snrTarget = max(guess.signalToNoise, snrPrevious)

            // This check is required to prevent low state go too close to 0, causing too broad peaks
            if (lowState.mean < guess.lowMin) {
                LOG.info("Low state mean ${lowState.mean} < ${guess.lowMin}, fixing...")
                outOfLowerNoise = true
                lowState.mean = guess.lowMin
                lowState.failures = estimateFailuresUsingMoments(
                    lowState.mean,
                    max(lowState.mean * SPAN_HMM_NB_VAR_MEAN_MULTIPLIER, lowState.variance)
                )
                updatedLow = true
            }

            // This check is required mostly for narrow marks to guard decent signal-to-noise ratio
            if (snrPrevious < snrTarget || updatedLow) {
                if (snrPrevious < snrTarget) {
                    LOG.info("Signal-to-noise ratio $snrPrevious < ${snrTarget}, fixing...")
                    outOfSignalToNoiseRatioRangeDown = true
                } else {
                    LOG.info("Updating high state mean, snr = $snrPrevious...")
                }
                highState.mean = max(highState.mean, lowState.mean * guess.signalToNoise)
                highState.failures = estimateFailuresUsingMoments(
                    highState.mean,
                    max(highState.mean * SPAN_HMM_NB_VAR_MEAN_MULTIPLIER, highState.variance)
                )
                updatedHigh = true
            }

            if (updatedLow) {
                lowState.updateTransients()
            }
            if (updatedHigh) {
                highState.updateTransients()
            }
        }
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 3

        // Will be updated during guess step
        lateinit var guess: NegBinUtil.Guess

        private val LOG: Logger = LoggerFactory.getLogger(NB2ZHMM::class.java)

        fun fitter(hmmEstimateSNR: Double, hmmLow: Double) = object : Fitter<NB2ZHMM> {
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
                guess = guessByData(positiveCoverage(preprocessed), 2,
                    estimateSNRFraction = hmmEstimateSNR,
                    estimateLowMinThreshold = hmmLow
                )
                return NB2ZHMM(guess.means, guess.failures)
            }
        }
    }
}
