package org.jetbrains.bio.statistics.hmm

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

        fun fitter() = object : Fitter<MLFreeNBHMM> {
            override fun guess(preprocessed: Preprocessed<DataFrame>, title: String,
                               threshold: Double, maxIter: Int): MLFreeNBHMM =
                    guess(listOf(preprocessed), title, threshold, maxIter)

            override fun guess(preprocessed: List<Preprocessed<DataFrame>>, title: String,
                               threshold: Double, maxIter: Int): MLFreeNBHMM {
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.toIntArray()
                val mean = emissions.average()
                check(mean != 0.0) { "Model can't be trained on empty coverage, exiting." }
                val sd = emissions.standardDeviation()
                val failures = NegativeBinomialDistribution
                        .estimateFailuresUsingMoments(mean, sd * sd)
                return MLFreeNBHMM(mean * .5, mean * 2.0, failures)
            }
        }
    }
}