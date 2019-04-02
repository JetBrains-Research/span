package org.jetbrains.bio.statistics.hmm

import org.apache.log4j.Logger
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.viktor.F64Array

/**
 * A zero-inflated HMM with univariate Negative Binomial emissions.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 25/05/15
 */
class MLFreeNBHMM(meanLow: Double, meanHigh: Double, failures: Double) : MLFreeHMM(3, 1) {
    private val zeroEmission: ConstantIntegerEmissionScheme = ConstantIntegerEmissionScheme(0)
    private val negBinEmissionSchemes: Array<NegBinEmissionScheme> = arrayOf(
            NegBinEmissionScheme(meanLow, failures),
            NegBinEmissionScheme(meanHigh, failures))

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

    val means: F64Array get() = F64Array(2) { negBinEmissionSchemes[it].mean }

    val failures: F64Array get() = F64Array(2) { negBinEmissionSchemes[it].failures }

    val successProbabilities: F64Array get() = F64Array(2) { negBinEmissionSchemes[it].successProbability }

    override fun toString(): String = toStringHelper()
            .add("means", means)
            .add("failures", failures)
            .toString()

    companion object {
        @Transient
        @JvmField
        val VERSION: Int = 1

        private val LOG = Logger.getLogger(MLFreeNBHMM::class.java)

        /**
         * This value is used to propose initial states for mean values of NOISE and SIGNAL states in the model.
         * Real life signal-to-noise ratio are generally in range 10-30.
         */
        const val SIGNAL_TO_NOISE_RATIO = 100.0

        fun fitter() = object : Fitter<MLFreeNBHMM> {
            override fun guess(preprocessed: Preprocessed<DataFrame>, title: String,
                               threshold: Double, maxIter: Int): MLFreeNBHMM =
                    guess(listOf(preprocessed), title, threshold, maxIter)

            override fun guess(preprocessed: List<Preprocessed<DataFrame>>, title: String,
                               threshold: Double, maxIter: Int): MLFreeNBHMM {
                // Filter out 0s, since they are covered by dedicated ZERO state
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.filter { it != 0 }.toIntArray()
                check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
                val mean = emissions.average()
                val sd = emissions.standardDeviation()
                val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, sd * sd)
                val meanLow = mean / Math.sqrt(SIGNAL_TO_NOISE_RATIO)
                val meanHigh = mean * Math.sqrt(SIGNAL_TO_NOISE_RATIO)
                LOG.debug("Guess emissions mean $mean\tsd $sd")
                LOG.debug("Guess init meanLow $meanLow\tmeanHigh $meanHigh\tfailures $fs")
                return MLFreeNBHMM(meanLow, meanHigh, fs)
            }
        }
    }
}
