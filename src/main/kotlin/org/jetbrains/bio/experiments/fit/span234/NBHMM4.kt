package org.jetbrains.bio.experiments.fit.span234

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLFreeHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import kotlin.math.max
import kotlin.math.pow

class NBHMM4(means: DoubleArray, failures: Double) : MLFreeHMM(4, 1) {
    init {
        check(means.size == 3) { "Expected 3 states" }
    }

    private val zeroEmission: ConstantIntegerEmissionScheme = ConstantIntegerEmissionScheme(0)
    private val negBinEmissionSchemes: Array<NegBinEmissionScheme> =
        Array(3) { NegBinEmissionScheme(means[it], failures) }

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
        NMHMMNZ.flipStatesIfNecessary(negBinEmissionSchemes, logPriorProbabilities, logTransitionProbabilities)
    }

    val means: F64Array get() = F64Array(2) { negBinEmissionSchemes[it].mean }

    val failures: F64Array get() = F64Array(2) { negBinEmissionSchemes[it].failures }

    val successProbabilities: F64Array get() = F64Array(2) { negBinEmissionSchemes[it].successProbability }

    override fun toString(): String = toStringHelper()
        .add("means", means)
        .add("failures", failures)
        .toString()

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 1

        private val LOG = LoggerFactory.getLogger(NBHMM4::class.java)

        fun fitter() = object : Fitter<NBHMM4> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM4 =
                guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>, title: String,
                threshold: Double, maxIter: Int, attempt: Int
            ): NBHMM4 {
                // Filter out 0s, since they are covered by dedicated ZERO state
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.filter { it != 0 }.toIntArray()
                check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
                val mean = emissions.average()
                val sd = emissions.standardDeviation()
                val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, sd * sd)
                val snr = MLFreeNBHMM.signalToNoise(attempt)
                val means = DoubleArray(3) {mean / snr.pow((3 / 2 - it).toDouble()) }
                LOG.debug("Guess $attempt emissions $means, failures $fs")
                return NBHMM4(means, fs)
            }
        }

        /**
         *
         * This value is used to propose initial states for mean values of LOW and HIGH states in the model.
         * Good experiment signal-to-noise ratio is generally 10-30.
         *
         * Use multiplier sequence depending on attempt number 1, 2, 1/2, 4, 1/4, etc.
         * Used for multistart
         */
        fun signalToNoise(attempt: Int, snr: Double = 20.0) =
            max(1.1, snr * 2.0.pow(((attempt + 1) / 2) * (if (attempt % 2 == 1) 1.0 else -1.0)))
    }
}
