package org.jetbrains.bio.span.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.statistics.regression.PoissonRegressionEmissionScheme
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.EmissionScheme
import org.jetbrains.bio.statistics.mixture.MLFreeMixture
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.asF64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import kotlin.math.exp

/**
 * Mixture of 3 components:
 * 0 - zero-inflated component
 * 1 - LOW (poisson with small parameters)
 * 2 - HIGH (poisson with high parameters)
 *
 * @author Elena Kartysheva
 * @date 5/25/19
 */
class PoissonRegression2Mixture(
    weights: F64Array,
    covariateLabels: List<String>,
    regressionCoefficients: Array<DoubleArray>
) : MLFreeMixture(numComponents = 3, numDimensions = 1, weights = weights) {

    private val zeroEmission = ConstantIntegerEmissionScheme(0)
    private val regressionEmissionSchemes = arrayOf(
        PoissonRegressionEmissionScheme(
            covariateLabels = covariateLabels,
            regressionCoefficients = regressionCoefficients[0]
        ),
        PoissonRegressionEmissionScheme(
            covariateLabels = covariateLabels,
            regressionCoefficients = regressionCoefficients[1]
        )
    )

    val signalToNoise
        get() = exp(
            regressionEmissionSchemes[1].regressionCoefficients[0] -
                    regressionEmissionSchemes[0].regressionCoefficients[0]
        )

    operator fun get(i: Int) = if (i == 0) zeroEmission else regressionEmissionSchemes[i - 1]

    operator fun set(i: Int, e: PoissonRegressionEmissionScheme) {
        require(i > 0)
        regressionEmissionSchemes[i - 1] = e
    }

    override fun getEmissionScheme(i: Int, d: Int): EmissionScheme {
        require(d == 0) { "Invalid dimension $d" }
        return get(i)
    }

    /**
     * We assume that the response vector is the integer-valued column 0,
     * and that the remaining columns include all the covariate labels as the double-valued covariates.
     */
    override fun fit(preprocessed: List<Preprocessed<DataFrame>>, title: String, threshold: Double, maxIterations: Int) {
        val data = DataFrame.rowBind(preprocessed.map { it.get() }.toTypedArray())
        super.fit(Preprocessed.of(data), title, threshold, maxIterations)
        flipStatesIfNecessary()
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        var VERSION = 2

        val LOG: Logger = LoggerFactory.getLogger(PoissonRegression2Mixture::class.java)

        fun fitter() = object : Fitter<PoissonRegression2Mixture> {
            /**
             * We assume that the response vector is the integer-valued column 0,
             * and that the remaining columns are the double-valued covariates.
             */
            override fun guess(
                preprocessed: Preprocessed<DataFrame>,
                title: String,
                threshold: Double,
                maxIterations: Int,
                attempt: Int
            ) = guess(listOf(preprocessed), title, threshold, maxIterations, attempt)

            /**
             * We assume that the response vector is the integer-valued column 0,
             * and that the remaining columns are the double-valued covariates.
             */
            override fun guess(
                preprocessed: List<Preprocessed<DataFrame>>,
                title: String,
                threshold: Double,
                maxIterations: Int,
                attempt: Int
            ): PoissonRegression2Mixture {
                require(attempt <= 1) { "Multistart is not supported." }
                // Filter out 0s, since they are covered by dedicated ZERO state
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.filter { it != 0 }.toIntArray()
                check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
                val df = preprocessed[0].get()
                return PoissonRegression2Mixture(
                    doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
                    df.labels.drop(1),
                    arrayOf(
                        DoubleArray(df.columnsNumber) { 0.0 },
                        DoubleArray(df.columnsNumber) { if (it == 0) 1.0 else 0.0 }
                    )
                )
            }
        }

        /**
         * Flip states in case when states with HIGH get lower mean than LOW
         */
        internal fun PoissonRegression2Mixture.flipStatesIfNecessary() {
            val lowScheme = this[1] as PoissonRegressionEmissionScheme
            val highScheme = this[2] as PoissonRegressionEmissionScheme
            if (weights[1] < weights[2]) {
                LOG.warn("After fitting the model, the weight of LOW state is lower than that of HIGH state.")
                LOG.warn("This usually indicates that the states were flipped during fitting. We will now flip them back.")
                this[2] = lowScheme
                this[1] = highScheme
                val tmp = logWeights[1]
                logWeights[1] = logWeights[2]
                logWeights[2] = tmp
            }
        }

    }
}
