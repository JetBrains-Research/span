package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array

/**
 * A zero-inflated HMM with univariate Negative Binomial emissions.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 25/05/15
 */
class NB2ZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) :
    FreeNBZHMM(nbMeans, nbFailures,
        priorProbabilities= NB2ZHMM_PRIORS,
        transitionProbabilities = NB2ZHMM_TRANSITIONS
    ) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 2

        // Proper initialisation makes convergence faster, more accurate
        val NB2ZHMM_PRIORS = F64Array.of(0.5, 0.45, 0.05)

        private val NB2ZHMM_TRANSITIONS_ = listOf(
            doubleArrayOf(0.99, 0.009, 0.001),
            doubleArrayOf(0.8, 0.19, 0.01),
            doubleArrayOf(0.01, 0.04, 0.95))

        val NB2ZHMM_TRANSITIONS = F64Array.invoke(3, 3) {i, j -> NB2ZHMM_TRANSITIONS_[i][j]}


        fun fitter() = object : Fitter<NB2ZHMM> {
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB2ZHMM = guess(listOf(preprocessed), title, threshold, maxIterations)

            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int
            ): NB2ZHMM {
                val (means, failures) = guess(preprocessed, 2)
                return NB2ZHMM(means, failures)
            }
        }
    }
}
