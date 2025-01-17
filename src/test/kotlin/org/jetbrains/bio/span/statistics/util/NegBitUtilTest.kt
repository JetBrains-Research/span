package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.distribution.Sampling
import org.junit.Before
import org.junit.Test
import kotlin.test.assertEquals

class NegBitUtilTest {

    @Before fun setUp() {
        Sampling.RANDOM_DATA_GENERATOR.reSeed(100)
    }

    @Test
    fun testGuessByData() {
        val sLow = NegBinEmissionScheme(1.0, 100.0).sampler()
        val sHigh = NegBinEmissionScheme(40.0, 1.0).sampler()

        val emissions = IntArray(1000) { sLow.asInt + 1 } + IntArray(100) { sHigh.asInt + 1 }
        emissions.shuffle()
        val (means, failures, lowMin, snrMin) = NegBinUtil.guessByData(
            emissions, 2
        )
        assertEquals(1.0, means[0], 1.0)
        assertEquals(34.0, means[1], 1.0)
        assertEquals(0.1, failures[0], .1)
        assertEquals(0.1, failures[1], .1)
        assertEquals(0.4, lowMin, .1)
        assertEquals(26.1, snrMin, 1.0)
    }

    @Test
    fun testGuessByData3() {
        val sLow = NegBinEmissionScheme(1.0, 100.0).sampler()
        val sMed = NegBinEmissionScheme(10.0, 30.0).sampler()
        val sHigh = NegBinEmissionScheme(40.0, 1.0).sampler()

        val emissions = IntArray(1000) { sLow.asInt + 1 } +
                IntArray(1000) { sMed.asInt + 1 } +
                IntArray(1000) { sHigh.asInt + 1 }
        emissions.shuffle()
        val (mean, _, _, _) = NegBinUtil.guessByData(
            emissions, 3
        )
        assertEquals(3.0, mean[0], 1.0)
        assertEquals(16.0, mean[1], 1.0)
        assertEquals(86.0, mean[2], 1.0)
    }

}