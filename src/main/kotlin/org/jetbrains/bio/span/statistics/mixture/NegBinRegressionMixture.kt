package org.jetbrains.bio.span.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.flipStatesIfNecessary
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.mixture.MLFreeMixture
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.span.statistics.regression.NegBinRegressionEmissionScheme
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import kotlin.math.exp

/**
 * Mixture of three components:
 * 0 - NULL (zero component)
 * 1 - LOW (negative binomial regression with low mean)
 * 2 - HIGH (negative binomial regression with high mean)
 *
 * @author Elena Kartysheva
 * @date 1/26/20
 */
class NegBinRegressionMixture(
    weights: F64Array,
    covariateLabels: List<String>,
    regressionCoefficients: Array<DoubleArray>
) : MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights) {

    private val zeroEmission = ConstantIntegerEmissionScheme(0)
    private val regressionEmissionSchemes = arrayOf(
        NegBinRegressionEmissionScheme(
            covariateLabels = covariateLabels,
            regressionCoefficients = regressionCoefficients[0],
            failures = 0.0
        ),
        NegBinRegressionEmissionScheme(
            covariateLabels = covariateLabels,
            regressionCoefficients = regressionCoefficients[1],
            failures = 0.0
        )
    )

    val signalToNoise
        get() = exp(
            regressionEmissionSchemes[1].regressionCoefficients[0] -
                    regressionEmissionSchemes[0].regressionCoefficients[0]
        )

    operator fun get(i: Int) = if (i == 0) zeroEmission else regressionEmissionSchemes[i - 1]

    operator fun set(i: Int, e: NegBinRegressionEmissionScheme) {
        require(i > 0)
        regressionEmissionSchemes[i - 1] = e
    }

    override fun getEmissionScheme(i: Int, d: Int): EmissionScheme {
        require(d == 0) { "Invalid dimension $d" }
        return get(i)
    }

    /**
     * We assume that the response vector is the integer-valued column 0,
     * and that the remaining columns include all of the covariate labels as the double-valued covariates.
     */
    override fun fit(preprocessed: List<Preprocessed<DataFrame>>, title: String, threshold: Double, maxIter: Int) {
        val data = DataFrame.rowBind(preprocessed.map { it.get() }.toTypedArray())
        super.fit(Preprocessed.of(data), title, threshold, maxIter)
        flipStatesIfNecessary()
    }

    companion object {
        @Transient
        @JvmField
        var VERSION = 1

        fun fitter() = object : Fitter<NegBinRegressionMixture> {
            /**
             * We assume that the response vector is the integer-valued column 0,
             * and that the remaining columns are the double-valued covariates.
             */
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIter: Int,
                attempt: Int
            ) = guess(listOf(preprocessed), title, threshold, maxIter, attempt)

            /**
             * We assume that the response vector is the integer-valued column 0,
             * and that the remaining columns are the double-valued covariates.
             */
            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIter: Int,
                attempt: Int
            ): NegBinRegressionMixture {
                // Filter out 0s, since they are covered by dedicated ZERO state
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.filter { it != 0 }.toIntArray()
                check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
                val df = preprocessed[0].get()
                return NegBinRegressionMixture(
                    doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
                    df.labels.drop(1),
                    arrayOf(
                        DoubleArray(df.columnsNumber) { 0.0 },
                        DoubleArray(df.columnsNumber) { if (it == 0) 1.0 else 0.0 }
                    )
                )
            }
        }
    }
}