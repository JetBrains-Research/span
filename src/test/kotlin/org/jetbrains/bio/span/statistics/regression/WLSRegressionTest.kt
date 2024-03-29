package org.jetbrains.bio.span.statistics.regression

import org.jetbrains.bio.Tests
import org.jetbrains.bio.viktor.asF64Array
import org.junit.Assert.assertEquals
import org.junit.Test
import kotlin.test.assertFails
import kotlin.test.assertTrue

class WLSRegressionTest {

    private val matrix123456 = arrayOf(
        doubleArrayOf(1.0, 2.0, 3.0),
        doubleArrayOf(4.0, 5.0, 6.0)
    )

    @Test
    fun testDesignMatrix() {
        // column size mismatch
        assertFails {
            WLSRegression.designMatrix(arrayOf(doubleArrayOf(1.0), doubleArrayOf(2.0, 3.0)))
        }
        // empty matrix
        assertFails {
            WLSRegression.designMatrix(arrayOf())
        }
        assertTrue(
            arrayOf(
                doubleArrayOf(1.0, 1.0, 1.0),
                doubleArrayOf(1.0, 2.0, 3.0),
                doubleArrayOf(4.0, 5.0, 6.0)
            ).contentDeepEquals(WLSRegression.designMatrix(matrix123456)), "Array contents differ"
        )
    }

    @Test
    fun testCalculateEta() {
        assertFails {
            WLSRegression.calculateEta(emptyArray(), doubleArrayOf(1.0))
        }
        assertFails {
            WLSRegression.calculateEta(matrix123456, doubleArrayOf())
        }
        assertFails {
            WLSRegression.calculateEta(matrix123456, doubleArrayOf(1.0, 2.0, 3.0))
        }
        assertEquals(
            doubleArrayOf(39.0, 54.0, 69.0).asF64Array(),
            WLSRegression.calculateEta(matrix123456, doubleArrayOf(7.0, 8.0))
        )
    }

    @Test
    fun testCalculateBeta() {
        assertFails {
            WLSRegression.calculateBeta(emptyArray(), doubleArrayOf(1.0).asF64Array(), doubleArrayOf(2.0).asF64Array())
        }
        assertFails {
            WLSRegression.calculateBeta(matrix123456, doubleArrayOf().asF64Array(), doubleArrayOf().asF64Array())
        }
        assertFails {
            WLSRegression.calculateBeta(
                matrix123456,
                doubleArrayOf(1.0, 2.0).asF64Array(),
                doubleArrayOf(3.0, 4.0, 5.0).asF64Array()
            )
        }
        assertFails {
            WLSRegression.calculateBeta(
                matrix123456,
                doubleArrayOf(1.0, 2.0, 3.0).asF64Array(),
                doubleArrayOf(4.0, 5.0).asF64Array()
            )
        }
        // data and reference values generated by R:
        // > X = data.frame(x1=runif(4), x2=runif(4), y=runif(4))
        // > w = runif(4)
        // > lm(y ~ x1 + x2, data=X, weights=w)
        val x = arrayOf(
            doubleArrayOf(0.0130203, 0.5831732, 0.1911626, 0.7604341),
            doubleArrayOf(0.76642767, 0.03048108, 0.52236819, 0.65585551)
        )
        val y = doubleArrayOf(0.5986931, 0.8171509, 0.5662360, 0.4927709).asF64Array()
        val weights = doubleArrayOf(0.180464329, 0.825891734, 0.008628437, 0.696381745).asF64Array()
        val expectedBeta = doubleArrayOf(0.95045, -0.20516, -0.45993)
        val actualBeta = WLSRegression.calculateBeta(WLSRegression.designMatrix(x), y, weights)
        Tests.assertEquals(expectedBeta, actualBeta, 1E-5)
    }
}
