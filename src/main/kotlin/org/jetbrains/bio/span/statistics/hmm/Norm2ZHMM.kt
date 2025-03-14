package org.jetbrains.bio.span.statistics.hmm

import org.apache.commons.logging.LogFactory
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_ESTIMATE_LOW
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_HMM_MAX_MEAN_TO_STD
import org.jetbrains.bio.span.statistics.emission.NormalEmissionScheme
import org.jetbrains.bio.span.statistics.hmm.FreeNBZHMM.Companion.positiveCoverage
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLFreeHMM
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.statistics.stochastic
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.max
import kotlin.math.pow

/**
 * Hidden Markov model with multidimensional integer-valued emissions.
 * Consists of ZERO state and 2 Gaussian Distributions
 */
class Norm2ZHMM(
    means: DoubleArray, variances: DoubleArray,
    priorProbabilities: F64Array = F64Array.stochastic(means.size + 1),
    transitionProbabilities: F64Array = F64Array.stochastic(means.size + 1, means.size + 1)
) : MLFreeHMM(means.size + 1, 1, priorProbabilities, transitionProbabilities) {

    private val zeroEmission: ConstantIntegerEmissionScheme = ConstantIntegerEmissionScheme(0)
    private val normalEmissionSchemes: Array<NormalEmissionScheme> =
        Array(means.size) { NormalEmissionScheme(means[it], variances[it]) }

    override fun getEmissionScheme(i: Int, d: Int): IntegerEmissionScheme {
        return if (i == 0) zeroEmission else normalEmissionSchemes[i - 1]
    }

    operator fun get(i: Int): IntegerEmissionScheme {
        return if (i == 0) zeroEmission else normalEmissionSchemes[i - 1]
    }

    operator fun set(i: Int, e: NormalEmissionScheme) {
        if (i == 0) {
            throw IllegalArgumentException()
        } else {
            normalEmissionSchemes[i - 1] = e
        }
    }

    override fun updateParameters(df: DataFrame, gammas: F64Array) {
        super.updateParameters(df, gammas)

        // You wanna keep em separated!
        for (d in 0 until numDimensions) {
            val lowState = getEmissionScheme(1, d) as NormalEmissionScheme
            val highState = getEmissionScheme(2, d) as NormalEmissionScheme

            val snrPrevious = highState.mean / lowState.mean
            var updated = false

            // This check is required to prevent low state go too close to 0, causing too broad peaks
            if (lowState.mean < guess.lowMin) {
                LOG.warn("Low state mean ${lowState.mean} < ${guess.lowMin}, fixing...")
                lowState.mean = guess.lowMin
                lowState.variance = max(lowState.variance, lowState.mean)
                updated = true
            }

            val snr = highState.mean / lowState.mean
            val snrTarget = max(guess.signalToNoise, snrPrevious)

            // This check is required mostly for narrow marks to guard decent signal-to-noise ratio
            if (snr < snrTarget) {
                LOG.warn("Signal-to-noise ratio $snr < ${snrTarget}, fixing...")
                highState.mean = lowState.mean * snrTarget
                updated = true
            }

            if (updated) {
                lowState.updateTransients()
                highState.updateTransients()
            }
        }
    }

    val means: F64Array get() = F64Array(normalEmissionSchemes.size) { normalEmissionSchemes[it].mean }

    val variances: F64Array get() = F64Array(normalEmissionSchemes.size) { normalEmissionSchemes[it].variance }

    override fun toString(): String = toStringHelper()
        .add("means", means)
        .add("variances", variances)
        .toString()

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1

        val LOG = LogFactory.getLog(Norm2ZHMM::class.java)

        @Suppress("ArrayInDataClass")
        data class Guess(
            val means: DoubleArray,
            val variances: DoubleArray,
            val lowMin: Double,
            val signalToNoise: Double
        )

        // Will be updated during guess step
        lateinit var guess: Guess


        fun guess(preprocessed: List<Preprocessed<DataFrame>>, n: Int): Guess {
            val emissions = preprocessed.flatMap {
                it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
            }.toIntArray()
            return guessByData(emissions, n)
        }


        /**
         * Guess initial approximation for mixture / HMM given data for several Normal Distributions
         */
        private fun guessByData(
            emissions: IntArray,
            n: Int,
            estimateSNRFraction: Double = SPAN_DEFAULT_HMM_ESTIMATE_SNR,
            estimateLowFraction: Double = SPAN_HMM_ESTIMATE_LOW,
            estimateLowMinThreshold: Double = SPAN_DEFAULT_HMM_LOW_THRESHOLD,
            maxMeanToStd: Double = SPAN_HMM_MAX_MEAN_TO_STD
        ): Guess {
            require(estimateSNRFraction + estimateLowFraction < 1.0) {
                "Estimate SNR fraction is too high $estimateSNRFraction"
            }
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

            // Most likely close to 1, but in case when all the coverage is extremely high this helps
            val lowEmissions = IntArray((emissions.size * estimateLowFraction).toInt()) {
                emissions[emissions.size - it - 1]
            }
            val meanL = lowEmissions.average()
            val sdL = lowEmissions.standardDeviation()
            LOG.debug("Low $estimateLowFraction emissions mean $meanL\t std $sdL")

            val signalToNoiseRatio = meanH / meanL

            // Make means proportional to snr, minimal should be equal to the noise state
            val means = DoubleArray(n) { meanL * signalToNoiseRatio.pow(it.toDouble() / (n - 1)) }

            val snrMin = if (estimateSNRFraction > 0) signalToNoiseRatio else 0.0
            LOG.debug("Minimal signal to noise ratio: $snrMin")
            val lowMin = meanL * estimateLowMinThreshold
            LOG.debug("Minimal low state: $lowMin")

            return Guess(
                means = means,
                variances = doubleArrayOf(sd * sd, sdH * sdH),
                lowMin = lowMin,
                signalToNoise = snrMin
            )
        }

        fun fitter(hmmEstimateSNR: Double, hmmLow: Double) = object : Fitter<Norm2ZHMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): Norm2ZHMM = guess(listOf(preprocessed), title, threshold, maxIterations)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): Norm2ZHMM {
                guess = guessByData(
                    positiveCoverage(preprocessed), 2,
                    estimateSNRFraction = hmmEstimateSNR,
                    estimateLowMinThreshold = hmmLow
                )
                return Norm2ZHMM(guess.means, guess.variances)
            }
        }


    }
}

