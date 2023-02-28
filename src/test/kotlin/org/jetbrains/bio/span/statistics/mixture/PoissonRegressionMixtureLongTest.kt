package org.jetbrains.bio.span.statistics.mixture

import org.jetbrains.bio.Tests
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.span.statistics.regression.PoissonRegressionEmissionScheme
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.random.Random


class PoissonRegressionMixtureLongTest {
    @Test
    fun sampleAndFit() {
        val original = PoissonRegression2Mixture(
            F64Array.of(0.2, 0.7, 0.1), listOf("foo", "bar"),
            arrayOf(
                doubleArrayOf(-1.0, 1.0, 1.0),
                doubleArrayOf(1.0, 1.0, 2.0)
            )
        )
        val rows = 500_000 // empirically found to produce sufficient re-learn precision
        val random = Random(1234)
        val df = DataFrame().with("y", IntArray(rows))
            .with("foo", DoubleArray(rows) { random.nextDouble() })
            .with("bar", DoubleArray(rows) { random.nextDouble() })
        original.sample(df, intArrayOf(0))
        val fitted = PoissonRegression2Mixture.fitter().fit(
            Preprocessed.of(df),
            title = "test",
            threshold = SpanPeakCallingExperiment.SPAN_FIT_THRESHOLD,
            maxIterations = SpanPeakCallingExperiment.SPAN_FIT_MAX_ITERATIONS
        )
        Tests.assertEquals(
            (original[1] as PoissonRegressionEmissionScheme).regressionCoefficients,
            (fitted[1] as PoissonRegressionEmissionScheme).regressionCoefficients,
            5E-2
        )
        Tests.assertEquals(
            (original[2] as PoissonRegressionEmissionScheme).regressionCoefficients,
            (fitted[2] as PoissonRegressionEmissionScheme).regressionCoefficients,
            5E-2
        )
        Tests.assertEquals(original.weights.toDoubleArray(), fitted.weights.toDoubleArray(), 5E-2)
    }
}