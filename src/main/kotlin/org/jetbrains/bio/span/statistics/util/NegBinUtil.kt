package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SIGNAL_TO_NOISE
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution.Companion.estimateFailuresUsingMoments
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
        signalToNoiseRatio: Double = SPAN_DEFAULT_SIGNAL_TO_NOISE
    ): Pair<DoubleArray, DoubleArray> {
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        val vars = sd * sd
        LOG.debug("Emissions mean $mean\t variance $vars")
        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val fs = estimateFailuresUsingMoments(mean, max(mean * 1.1, vars))
        LOG.debug("Failures {}", fs)
        // Make means proportional to snr.
        // Minimal should be equal to the genome-wide average, since we don't expect too much signal
        val means = DoubleArray(n) { mean * signalToNoiseRatio.pow(it.toDouble() / n) }
        val failures = DoubleArray(n) { fs }
        LOG.debug("Guess NegBinEmissionScheme means ${means.joinToString(",")}, failures ${failures.joinToString(",")}")
        return means to failures
    }
}