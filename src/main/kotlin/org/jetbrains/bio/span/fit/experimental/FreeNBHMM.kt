package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.span.statistics.util.NegBinUtil
import org.jetbrains.bio.span.statistics.util.NegBinUtil.guessByData
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLFreeHMM
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory

/**
 * Abstract hidden Markov model with multidimensional integer-valued emissions.
 * Consists of abstract number of Negative binomial states.
 */
open class FreeNBHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : MLFreeHMM(nbMeans.size, 1) {
    private val negBinEmissionSchemes: Array<NegBinEmissionScheme> =
        Array(nbMeans.size) { NegBinEmissionScheme(nbMeans[it], nbFailures[it]) }

    override fun getEmissionScheme(i: Int, d: Int): IntegerEmissionScheme {
        return negBinEmissionSchemes[i]
    }

    override fun fit(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String,
        threshold: Double,
        maxIterations: Int
    ) {
        super.fit(preprocessed, title, threshold, maxIterations)
        flipStatesIfNecessary(negBinEmissionSchemes, logPriorProbabilities, logTransitionProbabilities)
    }

    val means: F64Array get() = F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].mean }

    val failures: F64Array get() = F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].failures }

    val successProbabilities: F64Array
        get() =
            F64Array(negBinEmissionSchemes.size) { negBinEmissionSchemes[it].successProbability }

    override fun toString(): String = toStringHelper()
        .add("means", means)
        .add("failures", failures)
        .toString()

    companion object {
        private val LOG = LoggerFactory.getLogger(NB2HMM::class.java)

        fun guess(preprocessed: List<Preprocessed<DataFrame>>, n: Int): NegBinUtil.Guess {
            val emissions = preprocessed.flatMap {
                it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
            }.toIntArray()
            return guessByData(emissions, n)
        }

        /**
         * Flip states in case when states with HIGH get lower mean than LOW
         */
        fun flipStatesIfNecessary(
            negBinEmissionSchemes: Array<NegBinEmissionScheme>,
            logPriorProbabilities: F64Array,
            logTransitionProbabilities: F64Array
        ) {
            for (i in negBinEmissionSchemes.indices) {
                for (j in i + 1 until negBinEmissionSchemes.size) {
                    val lowScheme = negBinEmissionSchemes[i]
                    val highScheme = negBinEmissionSchemes[j]
                    val meanLow = lowScheme.mean
                    val meanHigh = highScheme.mean
                    val pLow = lowScheme.successProbability
                    val pHigh = highScheme.successProbability
                    val meanFlipped = meanLow > meanHigh
                    if (meanFlipped) {
                        LOG.debug(
                            "After fitting the model, mean emission in LOW state ($meanLow) is higher than " +
                                    "mean emission in HIGH state ($meanHigh)."
                        )
                    }
                    val pFlipped = pLow > pHigh
                    if (pFlipped) {
                        LOG.debug(
                            "After fitting the model, emission's parameter p in LOW state ($pLow) is higher than " +
                                    "emission's parameter p in HIGH state ($pHigh)."
                        )
                    }
                    if (meanFlipped && pFlipped) {
                        LOG.debug("Updated flipped states.")
                        negBinEmissionSchemes[i] = highScheme
                        negBinEmissionSchemes[j] = lowScheme
                        transitionProbabilityFlip(
                            i, j,
                            negBinEmissionSchemes, logPriorProbabilities, logTransitionProbabilities
                        )
                    } else if (meanFlipped || pFlipped) {
                        LOG.warn("Low quality of data detected after fitting the model.")
                    }
                }
            }
        }

        private fun transitionProbabilityFlip(
            state1: Int, state2: Int,
            negBinEmissionSchemes: Array<NegBinEmissionScheme>,
            logPriorProbabilities: F64Array,
            logTransitionProbabilities: F64Array
        ) {
            for (i in negBinEmissionSchemes.indices) {
                val tmp = logTransitionProbabilities[i, state1]
                logTransitionProbabilities[i, state1] = logTransitionProbabilities[i, state2]
                logTransitionProbabilities[i, state2] = tmp
            }
            for (j in negBinEmissionSchemes.indices) {
                val tmp = logTransitionProbabilities[state1, j]
                logTransitionProbabilities[state1, j] = logTransitionProbabilities[state2, j]
                logTransitionProbabilities[state2, j] = tmp
            }
            val tmp = logPriorProbabilities[state1]
            logPriorProbabilities[state1] = logPriorProbabilities[state2]
            logPriorProbabilities[state2] = tmp
        }

    }
}
