package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.experimental.FreeNBHMM
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLFreeHMM
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import kotlin.math.floor
import kotlin.math.max
import kotlin.math.pow

open class FreeNBZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : MLFreeHMM(nbMeans.size + 1, 1) {

    private val zeroEmission: ConstantIntegerEmissionScheme = ConstantIntegerEmissionScheme(0)
    private val negBinEmissionSchemes: Array<NegBinEmissionScheme> =
        Array(nbMeans.size) { NegBinEmissionScheme(nbMeans[it], nbFailures[it]) }

    override fun getEmissionScheme(i: Int, d: Int): IntegerEmissionScheme {
        return if (i == 0) zeroEmission else negBinEmissionSchemes[i - 1]
    }

    operator fun get(i: Int): IntegerEmissionScheme {
        return if (i == 0) zeroEmission else negBinEmissionSchemes[i - 1]
    }

    operator fun set(i: Int, e: NegBinEmissionScheme) {
        if (i == 0) {
            throw IllegalArgumentException()
        } else {
            negBinEmissionSchemes[i - 1] = e
        }
    }

    override fun fit(preprocessed: List<Preprocessed<DataFrame>>, title: String, threshold: Double, maxIter: Int) {
        super.fit(preprocessed, title, threshold, maxIter)
        flipStatesIfNecessary()
    }


    val means: F64Array get() = F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].mean }

    val failures: F64Array get() = F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].failures }

    val successProbabilities: F64Array
        get() = F64Array(negBinEmissionSchemes.size) {
            negBinEmissionSchemes[it].successProbability
        }

    override fun toString(): String = toStringHelper()
        .add("means", means)
        .add("failures", failures)
        .toString()

    fun flipStatesIfNecessary() {
        FreeNBHMM.flipStatesIfNecessary(negBinEmissionSchemes, logPriorProbabilities, logTransitionProbabilities)
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1

        private val LOG = LoggerFactory.getLogger(FreeNBZHMM::class.java)


        fun guess(preprocessed: List<Preprocessed<DataFrame>>, n: Int, attempt: Int): Pair<DoubleArray, DoubleArray> {
            // Filter out 0s, since they are covered by dedicated ZERO state
            val emissions = preprocessed.flatMap {
                it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
            }.filter { it != 0 }.toIntArray()
            check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
            val mean = emissions.average()
            val sd = emissions.standardDeviation()
            LOG.debug("Emissions mean $mean\t variance ${sd * sd}")
            // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
            // Otherwise, failures will be set to +Inf and won't be updated during EM steps.
            val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, max(1.1 * mean, sd * sd))
            val snr = multiStartSignalToNoise(attempt)
            // Make means and failures inverse proportional
            val means = DoubleArray(n) { mean / snr.pow((n / 2 - it).toDouble()) }
            val failures = DoubleArray(n) { fs * snr.pow((n / 2 - it).toDouble()) }
            LOG.debug("Guess $attempt emissions ${means.joinToString(",")}, failures ${failures.joinToString(",")}")
            return means to failures
        }

        /**
         * Propose initial signal-to-noise ration in multi-start runs.
         * Good experiment signal-to-noise ratio is generally 10-30.
         * [snr] and [multiplier] as used to yield the sequence of values
         * snr, snr / multiplier, snr * multiplier, snr / multiplier^2, snr * multiplier^2, ...
         */
        fun multiStartSignalToNoise(attempt: Int, snr: Double = 20.0, multiplier: Double = 2.0, min: Double = 1.1) =
            max(min, snr * multiplier.pow(floor((attempt + 1) / 2.0) * (-1.0).pow(attempt + 1)))

    }
}
