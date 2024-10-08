package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.fit.SpanConstants.SPAN_INITIAL_SCORES_HIGH
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_INITIAL_SCORES_LOW
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_NB_VAR_MEAN_MULTIPLIER
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_NOISE_MULTIPLIER
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
        val varsH: Double,
        val varsL: Double,
        val noiseMean: Double,
        val signalToNoise: Double
    )

    /**
     * Guess initial approximation for mixture / HMM given data for several Negative Binomial distributions.
     */
    fun guessByData(
        emissions: IntArray,
        n: Int,
    ): Guess {
        val mean = emissions.average()
        val sd = emissions.standardDeviation()
        val vars = sd * sd
        LOG.debug("All emissions mean $mean\t sd $sd")

        emissions.sortDescending()
        val highEmissions =
            IntArray((emissions.size * SPAN_INITIAL_SCORES_HIGH).toInt()) { emissions[it] }
        val meanH = highEmissions.average()
        val sdH = highEmissions.standardDeviation()
        LOG.debug("High $SPAN_INITIAL_SCORES_HIGH emissions mean $meanH\t std $sdH")

        val lowEmissions = IntArray((emissions.size * SPAN_INITIAL_SCORES_LOW).toInt()) {
            emissions[emissions.size - it - 1]
        }
        val meanL = lowEmissions.average()
        val sdL = lowEmissions.standardDeviation()
        LOG.debug("Low $SPAN_INITIAL_SCORES_LOW emissions mean $meanL\t std $sdL")

        val signalToNoiseRatio = (meanH + 1e-10) / (meanL + 1e-10)
        LOG.debug("Signal to noise ratio: $signalToNoiseRatio")
        // Make means proportional to snr.
        // Minimal should be equal to the genome-wide average, since we don't expect too much signal
        // Generally mean ~ meanL, but in synthetic examples, mean >> meanL, so we add compensation
        //
        // Assigning too high means for states leads to huge failures value, which can lead to the
        // situation, when high state successProbability is less than low state, which is considered
        // as a bad quality data marker and leads to incorrect "inverse" peak calling.
        // The example can be seen on simulation k36me3 by chips
        val means = DoubleArray(n) {
            sqrt(mean * meanL) * signalToNoiseRatio.pow(it.toDouble() / (n - 1))
        }

        // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
        // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
        val failures = DoubleArray(n) {
            estimateFailuresUsingMoments(means[it], max(means[it] * SPAN_NB_VAR_MEAN_MULTIPLIER, vars))
        }
        LOG.debug(
            "Guess NegBinEmissionScheme means ${means.joinToString(",")}, " +
                    "failures ${failures.joinToString(",")}"
        )
        return Guess(
            means,
            failures,
            sdH * sdH,
            sdL * sdL,
            meanL * SPAN_NOISE_MULTIPLIER,
            signalToNoiseRatio
        )
    }
}