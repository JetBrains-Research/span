package org.jetbrains.bio.span.fit

import org.jetbrains.bio.span.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.span.statistics.hmm.ConstrainedNBZHMM.Companion.flipStatesIfNecessary
import org.jetbrains.bio.span.statistics.hmm.FreeNBZHMM
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.span.statistics.mixture.PoissonRegressionMixture
import org.jetbrains.bio.span.statistics.mixture.PoissonRegressionMixture.Companion.flipStatesIfNecessary
import org.jetbrains.bio.span.statistics.regression.PoissonRegressionEmissionScheme
import org.jetbrains.bio.viktor.F64Array
import org.junit.Test
import kotlin.math.ln
import kotlin.test.assertEquals
import kotlin.test.assertFails

/**
 * @author Alexey Dievsky
 * @date 14.12.2015
 *
 * In the names of the tests:
 * - "F" means free model (just one replicate),
 * - "C" means constrained enrichment (>1 replicate),
 * - "D" means difference (1+ replicate vs 1+ replicate),
 * - "M" means Poisson regression mixture (just one replicate).
 */

class SpanFlipStatesTest {

    @Test
    fun testFirstStateFlipD() {
        val model = ConstrainedNBZHMM(
            ZLHID.constraintMap(3, 2),
            doubleArrayOf(12.0, 8.0, 10.0, 1.0, 4.0, 3.0, 0.5, 0.25, 0.4, 8.0),
            doubleArrayOf(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.4, 2.0)
        )
        /* there is no simple way to create a model with the specified
         * transition probabilities, so we set them manually */
        val logTransitionProbabilities = initTransitionProbsD(model)
        model.flipStatesIfNecessary(3, 2)
        val means = model.means
        val ps = model.successProbabilities
        assertFlipped(0, 3, means, ps)
        assertFlipped(1, 4, means, ps)
        assertFlipped(2, 5, means, ps)
        assertNotFlipped(6, 8, means, ps)
        assertNotFlipped(7, 9, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
        assertEquals(logTransitionProbabilities[3, 3], ln(.5))
        assertEquals(logTransitionProbabilities[4, 4], ln(.5))
        assertEquals(logTransitionProbabilities[0, 2], ln(.2))
    }

    @Test
    fun testSecondStateFlipD() {
        val model = ConstrainedNBZHMM(
            ZLHID.constraintMap(2, 3),
            doubleArrayOf(5.0, 0.25, 4.0, 8.0, 12.0, 8.0, 10.0, 1.0, 4.0, 3.0),
            doubleArrayOf(15.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        /* there is no simple way to create a model with the specified
         * transition probabilities, so we set them manually */
        val logTransitionProbabilities = initTransitionProbsD(model)
        model.flipStatesIfNecessary(2, 3)
        val means = model.means
        val ps = model.successProbabilities
        assertNotFlipped(0, 2, means, ps)
        assertNotFlipped(1, 3, means, ps)
        assertFlipped(4, 7, means, ps)
        assertFlipped(5, 8, means, ps)
        assertFlipped(6, 9, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
        assertEquals(logTransitionProbabilities[3, 3], ln(.5))
        assertEquals(logTransitionProbabilities[4, 4], ln(.5))
        assertEquals(logTransitionProbabilities[0, 1], ln(.1))
    }

    @Test
    fun testTwoStatesFlipD() {
        val model = ConstrainedNBZHMM(
            ZLHID.constraintMap(2, 3),
            doubleArrayOf(4.0, 8.0, 0.5, 0.25, 12.0, 8.0, 10.0, 1.0, 4.0, 3.0),
            doubleArrayOf(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        val logTransitionProbabilities = initTransitionProbsD(model)
        model.flipStatesIfNecessary(2, 3)
        val means = model.means
        val ps = model.successProbabilities
        assertFlipped(0, 2, means, ps)
        assertFlipped(1, 3, means, ps)
        assertFlipped(4, 7, means, ps)
        assertFlipped(5, 8, means, ps)
        assertFlipped(6, 9, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
        assertEquals(logTransitionProbabilities[3, 3], ln(.5))
        assertEquals(logTransitionProbabilities[4, 4], ln(.5))
    }

    private fun initTransitionProbsD(model: ConstrainedNBZHMM): F64Array {
        val logTransitionProbabilities = model.logTransitionProbabilities
        logTransitionProbabilities[0, 0] = ln(.5)
        logTransitionProbabilities[0, 1] = ln(.2)
        logTransitionProbabilities[0, 2] = ln(.1)
        logTransitionProbabilities[0, 3] = ln(.1)
        logTransitionProbabilities[0, 4] = ln(.1)
        logTransitionProbabilities[1, 0] = ln(.1)
        logTransitionProbabilities[1, 1] = ln(.5)
        logTransitionProbabilities[1, 2] = ln(.2)
        logTransitionProbabilities[1, 3] = ln(.1)
        logTransitionProbabilities[1, 4] = ln(.1)
        logTransitionProbabilities[2, 0] = ln(.2)
        logTransitionProbabilities[2, 1] = ln(.1)
        logTransitionProbabilities[2, 2] = ln(.5)
        logTransitionProbabilities[2, 3] = ln(.1)
        logTransitionProbabilities[2, 4] = ln(.1)
        logTransitionProbabilities[3, 0] = ln(.1)
        logTransitionProbabilities[3, 1] = ln(.1)
        logTransitionProbabilities[3, 2] = ln(.1)
        logTransitionProbabilities[3, 3] = ln(.5)
        logTransitionProbabilities[3, 4] = ln(.2)
        logTransitionProbabilities[4, 0] = ln(.2)
        logTransitionProbabilities[4, 1] = ln(.1)
        logTransitionProbabilities[4, 2] = ln(.1)
        logTransitionProbabilities[4, 3] = ln(.1)
        logTransitionProbabilities[4, 4] = ln(.5)
        return logTransitionProbabilities
    }

    @Test
    fun testNoStatesFlipD() {
        val model = ConstrainedNBZHMM(
            ZLHID.constraintMap(2, 3),
            doubleArrayOf(0.5, 8.0, 4.0, 8.0, 1.0, 4.0, 3.0, 12.0, 8.0, 10.0),
            doubleArrayOf(2.0, 12.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        val logTransitionProbabilities = initTransitionProbsD(model)
        model.flipStatesIfNecessary(2, 3)
        val means = model.means
        val ps = model.successProbabilities
        assertNotFlipped(0, 2, means, ps)
        assertNotFlipped(1, 3, means, ps)
        assertNotFlipped(4, 7, means, ps)
        assertNotFlipped(5, 8, means, ps)
        assertNotFlipped(6, 9, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
        assertEquals(logTransitionProbabilities[3, 3], ln(.5))
        assertEquals(logTransitionProbabilities[4, 4], ln(.5))
    }

    @Test
    fun testFittingErrorD() {
        val model = ConstrainedNBZHMM(
            ZLHID.constraintMap(2, 3),
            doubleArrayOf(4.0, 0.25, 0.5, 8.0, 1.0, 12.0, 3.0, 4.0, 8.0, 10.0),
            doubleArrayOf(2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        initTransitionProbsD(model)
        assertFails { model.flipStatesIfNecessary(2, 3) }
    }

    @Test
    fun testStateFlipC() {
        val model = ConstrainedNBZHMM(
            ZLH.constraintMap(3),
            doubleArrayOf(12.0, 8.0, 10.0, 1.0, 4.0, 3.0),
            doubleArrayOf(2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        val logTransitionProbabilities = initTransitionProbsC(model)
        model.flipStatesIfNecessary(3)
        val means = model.means
        val ps = model.successProbabilities
        assertFlipped(0, 3, means, ps)
        assertFlipped(1, 4, means, ps)
        assertFlipped(2, 5, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
    }

    private fun initTransitionProbsC(model: ConstrainedNBZHMM): F64Array {
        val logTransitionProbabilities = model.logTransitionProbabilities
        logTransitionProbabilities[0, 0] = ln(.5)
        logTransitionProbabilities[0, 1] = ln(.2)
        logTransitionProbabilities[0, 2] = ln(.3)
        logTransitionProbabilities[1, 0] = ln(.3)
        logTransitionProbabilities[1, 1] = ln(.5)
        logTransitionProbabilities[1, 2] = ln(.2)
        logTransitionProbabilities[2, 0] = ln(.2)
        logTransitionProbabilities[2, 1] = ln(.3)
        logTransitionProbabilities[2, 2] = ln(.5)
        return logTransitionProbabilities
    }

    @Test
    fun testNoStatesFlipC() {
        val model = ConstrainedNBZHMM(
            ZLH.constraintMap(3),
            doubleArrayOf(15.0, 4.0, 3.0, 12.0, 8.0, 10.0),
            doubleArrayOf(15.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        val logTransitionProbabilities = initTransitionProbsC(model)
        model.flipStatesIfNecessary(3)
        val means = model.means
        val ps = model.successProbabilities
        assertNotFlipped(0, 3, means, ps)
        assertNotFlipped(1, 4, means, ps)
        assertNotFlipped(2, 5, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
    }

    @Test
    fun testFittingErrorC() {
        val model = ConstrainedNBZHMM(
            ZLH.constraintMap(3),
            doubleArrayOf(1.0, 8.0, 3.0, 12.0, 4.0, 10.0),
            doubleArrayOf(2.0, 2.0, 2.0, 2.0, 2.0, 2.0)
        )
        initTransitionProbsC(model)
        assertFails { model.flipStatesIfNecessary(3) }
    }

    @Test
    fun testStateFlipF() {
        val model = NB2ZHMM(doubleArrayOf(8.0, 3.0), doubleArrayOf(2.0, 2.0))
        val logTransitionProbabilities = initTransitionProbsF(model)
        model.flipStatesIfNecessary()
        val means = model.means
        val ps = model.successProbabilities
        assertFlipped(0, 1, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
    }

    private fun initTransitionProbsF(model: FreeNBZHMM): F64Array {
        val logTransitionProbabilities = model.logTransitionProbabilities
        logTransitionProbabilities[0, 0] = ln(.5)
        logTransitionProbabilities[0, 1] = ln(.2)
        logTransitionProbabilities[0, 2] = ln(.3)
        logTransitionProbabilities[1, 0] = ln(.3)
        logTransitionProbabilities[1, 1] = ln(.5)
        logTransitionProbabilities[1, 2] = ln(.2)
        logTransitionProbabilities[2, 0] = ln(.2)
        logTransitionProbabilities[2, 1] = ln(.3)
        logTransitionProbabilities[2, 2] = ln(.5)
        return logTransitionProbabilities
    }

    @Test
    fun testNoStatesFlipF() {
        val model = NB2ZHMM(doubleArrayOf(8.0, 3.0), doubleArrayOf(2.0, 2.0))
        val logTransitionProbabilities = initTransitionProbsF(model)
        model.flipStatesIfNecessary()
        val means = model.means
        val ps = model.successProbabilities
        assertNotFlipped(0, 1, means, ps)
        assertEquals(logTransitionProbabilities[0, 0], ln(.5))
        assertEquals(logTransitionProbabilities[1, 1], ln(.5))
        assertEquals(logTransitionProbabilities[2, 2], ln(.5))
    }

    private fun assertFlipped(low: Int, high: Int, means: F64Array, ps: F64Array) {
        assert(means[low] < means[high] && ps[low] < ps[high]) { "Expected states to flip" }
    }

    private fun assertNotFlipped(low: Int, high: Int, means: F64Array, ps: F64Array) {
        assert(means[low] < means[high] || ps[low] < ps[high]) { "Expected states not to flip" }
    }

    @Test
    fun testNoStatesFlipM() {
        val regressionCoefficients = arrayOf(
            doubleArrayOf(-4.0, 1.0),
            doubleArrayOf(-5.0, 2.0)
        )
        val model = PoissonRegressionMixture(
            F64Array.of(0.5, 0.4, 0.1),
            listOf("first", "second"),
            regressionCoefficients
        )
        model.flipStatesIfNecessary()
        assertCorrect(model.weights)
        assertEquals(regressionCoefficients[0], (model[1] as PoissonRegressionEmissionScheme).regressionCoefficients)
        assertEquals(regressionCoefficients[1], (model[2] as PoissonRegressionEmissionScheme).regressionCoefficients)
    }

    @Test
    fun testFlipM() {
        val regressionCoefficients = arrayOf(
            doubleArrayOf(-5.0, 2.0),
            doubleArrayOf(-4.0, 1.0)
        )
        val model = PoissonRegressionMixture(
            F64Array.of(0.5, 0.1, 0.4),
            listOf("first", "second"),
            regressionCoefficients
        )
        model.flipStatesIfNecessary()
        assertCorrect(model.weights)
        assertEquals(regressionCoefficients[1], (model[1] as PoissonRegressionEmissionScheme).regressionCoefficients)
        assertEquals(regressionCoefficients[0], (model[2] as PoissonRegressionEmissionScheme).regressionCoefficients)
    }

    private fun assertCorrect(weights: F64Array) {
        assert(weights[1] >= weights[2]) { "LOW weight is less than HIGH weight" }
    }

}