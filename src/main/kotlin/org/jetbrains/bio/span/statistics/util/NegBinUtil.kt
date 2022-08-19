package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.standardDeviation
import org.slf4j.LoggerFactory
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
        snr: Double = 10.0
    ): Pair<DoubleArray, DoubleArray> {
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        LOG.debug("Emissions mean $mean\t variance ${sd * sd}")
        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, max(1.1 * mean, sd * sd))
        // Make means and failures inverse proportional to snr.
        val means = DoubleArray(n) { mean * snr.pow((it.toDouble() + 1) / n) }
        // NegativeBinomialDistribution variance greater than mean guard
        val failures = DoubleArray(n) { max(means[it] * 1.1, fs * snr.pow((it.toDouble() - n) / n)) }
        LOG.debug("Guess emissions ${means.joinToString(",")}, failures ${failures.joinToString(",")}")
        return means to failures
    }
}