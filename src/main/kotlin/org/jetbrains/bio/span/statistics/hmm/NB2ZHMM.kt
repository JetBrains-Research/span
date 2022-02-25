package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.model.Fitter

/**
 * A zero-inflated HMM with univariate Negative Binomial emissions.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 25/05/15
 */
class NB2ZHMM(nbMeans: DoubleArray, nbFailures: DoubleArray) : FreeNBZHMM(nbMeans, nbFailures) {

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 2


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
