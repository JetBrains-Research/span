package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_SIGNAL_ESTIMATE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_MAX_MEAN_TO_STD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_NB_VAR_MEAN_MULTIPLIER
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution.Companion.estimateFailuresUsingMoments
import org.jetbrains.bio.statistics.standardDeviation
import org.slf4j.LoggerFactory
import kotlin.math.max
import kotlin.math.pow
import kotlin.math.sqrt

object NegBinUtil {

    private val LOG = LoggerFactory.getLogger(NegBinUtil::class.java)

    @Suppress("ArrayInDataClass")
    data class Guess(
        val means: DoubleArray,
        val failures: DoubleArray,
        val lowMean: Double,
        val lowVar: Double,
        val signalToNoise: Double
    )

    /**
     * Guess initial approximation for mixture / HMM given data for several Negative Binomial distributions.
     */
    fun guessByData(
        emissions: IntArray,
        n: Int,
        estimateSignalFraction: Double = SPAN_HMM_SIGNAL_ESTIMATE,
        maxMeanToStd: Double = SPAN_HMM_MAX_MEAN_TO_STD,
        eps: Double = 1e-10
    ): Guess {
        require(estimateSignalFraction < 0.5) { "Estimate signal is too high $estimateSignalFraction" }
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        val vars = sd * sd
        LOG.debug("All emissions mean $mean\t std $sd")
        emissions.sortDescending()

        val highEmissions =
            IntArray((emissions.size * estimateSignalFraction).toInt()) { emissions[it] }
        var meanH = highEmissions.average()
        val sdH = highEmissions.standardDeviation()
        LOG.debug("High $estimateSignalFraction emissions mean $meanH\t std $sdH")
        if (meanH > maxMeanToStd * (sdH + eps)) {
            LOG.warn("High mean / std > $maxMeanToStd, adjusting...")
            meanH = (sdH + eps) * maxMeanToStd
            LOG.debug("Adjusted high $estimateSignalFraction emissions mean $meanH\t std $sdH")
        }
        val lowEstimateFraction = (1.0 - estimateSignalFraction) / 2
        val lowEmissions =
            IntArray((emissions.size * lowEstimateFraction).toInt()) { emissions[emissions.size - it - 1] }
        val meanL = lowEmissions.average()
        val sdL = lowEmissions.standardDeviation()
        LOG.debug("Low $lowEstimateFraction emissions mean $meanL\t std $sdL")

        val signalToNoiseRatio = (meanH + eps) / (meanL + eps)
        LOG.debug("Signal to noise ratio: $signalToNoiseRatio")
        // Make means proportional to snr.
        // Minimal should be equal to the genome-wide average, since we don't expect too much signal
        // Generally mean ~ meanL, but in synthetic examples, mean >> meanL, so we add compensation
        val means = DoubleArray(n) {
            sqrt(meanL * mean) * signalToNoiseRatio.pow(it.toDouble() / (n - 1))
        }

        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val failures = DoubleArray(n) {
            estimateFailuresUsingMoments(means[it], max(means[it] * SPAN_HMM_NB_VAR_MEAN_MULTIPLIER, vars))
        }
        LOG.debug(
            "Guess NegBinEmissionScheme means ${means.joinToString(",")}, " +
                    "failures ${failures.joinToString(",")}"
        )
        return Guess(
            means,
            failures,
            meanL,
            // Generally sdL ~ sd, but in synthetic or low coverage examples, sdL=0, so we add compensation
            if (sdL > 0) sqrt(sdL * sdL * vars) else vars,
            signalToNoiseRatio
        )
    }
}