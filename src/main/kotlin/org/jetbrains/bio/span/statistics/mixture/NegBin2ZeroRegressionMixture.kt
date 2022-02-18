package org.jetbrains.bio.span.statistics.mixture

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.statistics.regression.NegBinRegressionEmissionScheme
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
 * Mixture of three components:
 * 0 - NULL (zero component)
 * 1 - LOW (negative binomial regression with low mean)
 * 2 - HIGH (negative binomial regression with high mean)
 *
 * @author Elena Kartysheva
 * @date 1/26/20
 */
class NegBin2ZeroRegressionMixture(
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

        val LOG: Logger = LoggerFactory.getLogger(NegBin2ZeroRegressionMixture::class.java)

        fun fitter() = object : Fitter<NegBin2ZeroRegressionMixture> {
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
            ): NegBin2ZeroRegressionMixture {
                require(attempt <= 1) { "Multistart not supported." }
                // Filter out 0s, since they are covered by dedicated ZERO state
                val emissions = preprocessed.flatMap {
                    it.get().let { df -> df.sliceAsInt(df.labels.first()).toList() }
                }.filter { it != 0 }.toIntArray()
                check(emissions.isNotEmpty()) { "Model can't be trained on empty coverage, exiting." }
                val df = preprocessed[0].get()
                return NegBin2ZeroRegressionMixture(
                    doubleArrayOf(1 / 3.0, 1 / 3.0, 1 / 3.0).asF64Array(),
                    df.labels.drop(1),
                    arrayOf(
                        DoubleArray(df.columnsNumber) { 0.0 },
                        DoubleArray(df.columnsNumber) { if (it == 0) 1.0 else 0.0 }
                    )
                )
            }
        }

        internal fun NegBin2ZeroRegressionMixture.flipStatesIfNecessary() {
            val lowScheme = this[1] as NegBinRegressionEmissionScheme
            val highScheme = this[2] as NegBinRegressionEmissionScheme
            if (weights[1] < weights[2]) {
                LOG.warn("After fitting the model, the wight of LOW state is lower than that of HIGH state.")
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