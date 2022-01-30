package org.jetbrains.bio.span.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.span.statistics.hmm.FreeNBZHMM
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.mixture.MLFreeMixture
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory

/**
 * Mixture of three components:
 * 0 - NULL (zero component)
 * 1 - LOW (negative binomial regression with low mean)
 * 2 - HIGH (negative binomial regression with high mean)
 *
 * @author Oleg Shpynov
 * @date 1/29/22
 */
class NegBinMixture(
    nbMeans: DoubleArray, nbFailures: DoubleArray,
    weights: F64Array,
) : MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights) {

    private val zeroEmission: ConstantIntegerEmissionScheme = ConstantIntegerEmissionScheme(0)
    private val negBinEmissionSchemes: Array<NegBinEmissionScheme> =
        Array(nbMeans.size) { NegBinEmissionScheme(nbMeans[it], nbFailures[it]) }

    override fun getEmissionScheme(i: Int, d: Int): IntegerEmissionScheme {
        require(d == 0) { "Invalid dimension $d" }
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

    val means: F64Array get() = F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].mean }

    val failures: F64Array get() = F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].failures }

    override fun fit(preprocessed: List<Preprocessed<DataFrame>>, title: String, threshold: Double, maxIter: Int) {
        val data = DataFrame.rowBind(preprocessed.map { it.get() }.toTypedArray())
        super.fit(Preprocessed.of(data), title, threshold, maxIter)
        flipStatesIfNecessary()
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        var VERSION = 2

        val LOG: Logger = LoggerFactory.getLogger(NegBinMixture::class.java)

        fun fitter() = object : Fitter<NegBinMixture> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIter: Int,
                attempt: Int
            ) = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIter: Int,
                attempt: Int
            ): NegBinMixture {
                val (means, failures) = FreeNBZHMM.guess(preprocessed, 2, attempt)
                return NegBinMixture(
                    means, failures,
                    doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array()
                )
            }
        }

        internal fun NegBinMixture.flipStatesIfNecessary() {
            val lowScheme = negBinEmissionSchemes[0]
            val highScheme = negBinEmissionSchemes[1]
            val meanLow = lowScheme.mean
            val meanHigh = highScheme.mean
            val pLow = lowScheme.successProbability
            val pHigh = highScheme.successProbability
            val meanFlipped = meanLow > meanHigh
            if (meanFlipped) {
                LOG.warn(
                    "After fitting the model, mean emission in LOW state ($meanLow) is higher than " +
                            "mean emission in HIGH state ($meanHigh)."
                )
            }
            val pFlipped = pLow > pHigh
            if (pFlipped) {
                LOG.warn(
                    "After fitting the model, emission's parameter p in LOW state ($pLow) is higher than " +
                            "emission's parameter p in HIGH state ($pHigh)."
                )
            }
            if (meanFlipped && pFlipped) {
                LOG.warn("This usually indicates that the states were flipped during fitting. We will now flip them back.")
                negBinEmissionSchemes[0] = highScheme
                negBinEmissionSchemes[1] = lowScheme
                val tmp = weights[1]
                weights[1] = weights[2]
                weights[2] = tmp
            } else if (meanFlipped || pFlipped) {
                LOG.warn("This is generally harmless, but could indicate low quality of data.")
            }
        }
    }
}