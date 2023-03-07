package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.experimental.FreeNBHMM
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.span.statistics.util.NegBinUtil.guessByData
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLFreeHMM
import org.jetbrains.bio.viktor.F64Array

/**
 * Abstract hidden Markov model with multidimensional integer-valued emissions.
 * Consists of ZERO state and abstract number of Negative binomial states.
 */
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

    override fun fit(
        preprocessed: List<Preprocessed<DataFrame>>,
        title: String,
        threshold: Double,
        maxIterations: Int
    ) {
        super.fit(preprocessed, title, threshold, maxIterations)
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

        fun guess(preprocessed: List<Preprocessed<DataFrame>>, n: Int): Pair<DoubleArray, DoubleArray> {
            // Filter out 0s, since they are covered by dedicated ZERO state
            val emissions = preprocessed.flatMap {
                it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
            }.filter { it != 0 }.toIntArray()
            check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
            return guessByData(emissions, n)
        }

    }
}
