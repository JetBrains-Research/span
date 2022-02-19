package org.jetbrains.bio.span.statistics.emission

import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.standardDeviation
import org.slf4j.LoggerFactory
import kotlin.math.floor
import kotlin.math.max
import kotlin.math.pow

object NegBinUtil {

    private val LOG = LoggerFactory.getLogger(NegBinUtil::class.java)

    /**
     * Guess initial approximation for mixture / HMM given data for several Negative Binomial distributions.
     */
    fun guessByData(
        emissions: IntArray,
        n: Int,
        attempt: Int
    ): Pair<DoubleArray, DoubleArray> {
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        LOG.debug("Emissions mean $mean\t variance ${sd * sd}")
        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, max(1.1 * mean, sd * sd))
        val snr = multiStartSignalToNoise(attempt)
        // Make means and failures inverse proportional
        val means = DoubleArray(n) { mean * snr.pow(((it + 1) / n).toDouble()) }
        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val failures = DoubleArray(n) { max(means[it] * 1.1, fs * snr.pow(((it - n) / n).toDouble())) }
        LOG.debug("Guess $attempt emissions ${means.joinToString(",")}, failures ${failures.joinToString(",")}")
        return means to failures
    }

    /**
     * Propose initial signal-to-noise ration in multi-start runs.
     * [snr] and [multiplier] as used to yield the sequence of values
     * snr, snr / multiplier, snr * multiplier, snr / multiplier^2, snr * multiplier^2, ...
     */
    fun multiStartSignalToNoise(attempt: Int, snr: Double = 10.0, multiplier: Double = 2.0, min: Double = 1.1) =
        max(min, snr * multiplier.pow(floor((attempt + 1) / 2.0) * (-1.0).pow(attempt + 1)))

}