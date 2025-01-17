package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_ESTIMATE_LOW
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_MAX_MEAN_TO_STD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_NB_VAR_MEAN_MULTIPLIER
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution.Companion.estimateFailuresUsingMoments
import org.jetbrains.bio.statistics.standardDeviation
import org.slf4j.LoggerFactory
import kotlin.math.max
import kotlin.math.pow

object NegBinUtil {

    private val LOG = LoggerFactory.getLogger(NegBinUtil::class.java)

    @Suppress("ArrayInDataClass")
    data class Guess(
        val means: DoubleArray,
        val failures: DoubleArray,
        val lowMin: Double,
        val signalToNoise: Double
    )

    /**
     * Guess initial approximation for mixture / HMM given data for several Negative Binomial distributions.
     */
    fun guessByData(
        emissions: IntArray,
        n: Int,
        estimateSNRFraction: Double = SPAN_DEFAULT_HMM_ESTIMATE_SNR,
        estimateLowFraction: Double = SPAN_HMM_ESTIMATE_LOW,
        estimateLowMinThreshold: Double = SPAN_DEFAULT_HMM_LOW_THRESHOLD,
        maxMeanToStd: Double = SPAN_HMM_MAX_MEAN_TO_STD
    ): Guess {
        require(estimateSNRFraction + estimateLowFraction < 1.0) { "Estimate SNR fraction is too high $estimateSNRFraction" }
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        val vars = sd * sd
        LOG.debug("All emissions mean $mean\t std $sd")
        emissions.sortDescending()

        val fraction = if (estimateSNRFraction > 0) estimateSNRFraction else SPAN_DEFAULT_HMM_ESTIMATE_SNR
        val highEmissions = IntArray((emissions.size * fraction).toInt()) { emissions[it] }
        var meanH = highEmissions.average()
        val sdH = highEmissions.standardDeviation()
        LOG.debug("High $fraction emissions mean $meanH\t std $sdH")
        if (meanH > maxMeanToStd * (sdH + 1e-10)) {
            LOG.warn("High mean / std > $maxMeanToStd, adjusting...")
            meanH = (sdH + 1e-10) * maxMeanToStd
            LOG.debug("Adjusted high $fraction emissions mean $meanH\t std $sdH")
        }

        val lowEmissions = IntArray((emissions.size * estimateLowFraction).toInt()) {
            emissions[emissions.size - it - 1]
        }
        val meanL = lowEmissions.average()
        val sdL = lowEmissions.standardDeviation()
        LOG.debug("Low $estimateLowFraction emissions mean $meanL\t std $sdL")

        val signalToNoiseRatio = meanH / meanL

        // Make means proportional to snr, minimal should be equal to the noise state
        val means = DoubleArray(n) { meanL * signalToNoiseRatio.pow(it.toDouble() / (n - 1)) }
        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val failures = DoubleArray(n) {
            estimateFailuresUsingMoments(mean, max(mean * SPAN_HMM_NB_VAR_MEAN_MULTIPLIER, vars))
        }
        LOG.debug(
            "Guess NegBinEmissionScheme means ${means.joinToString(",")}, " +
                    "failures ${failures.joinToString(",")}"
        )

        val snrMin = if (estimateSNRFraction > 0) signalToNoiseRatio else 0.0
        LOG.debug("Minimal signal to noise ratio: $snrMin")
        val lowMin = meanL * estimateLowMinThreshold
        LOG.debug("Minimal low state: $lowMin")

        return Guess(
            means = means,
            failures = failures,
            lowMin = lowMin,
            signalToNoise = snrMin
        )
    }
}