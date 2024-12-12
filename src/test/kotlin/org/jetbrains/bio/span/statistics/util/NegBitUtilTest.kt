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
        val sLow = NegBinEmissionScheme(0.5, 100.0).sampler()
        val sHigh = NegBinEmissionScheme(40.0, 1.0).sampler()

        val emissions = IntArray(10000) { sLow.asInt + 1 } + IntArray(100) { sHigh.asInt + 1 }
        emissions.shuffle()
        val (means0, failures0, _, _) = NegBinUtil.guessByData(
            emissions, 2
        )
        assertEquals(1.3, means0[0], .1)
        assertEquals(8.7, means0[1], .1)
        assertEquals(0.1, failures0[0], .1)
        assertEquals(4.8, failures0[1], .1)
    }

    @Test
    fun testGuessByData3() {
        val sLow = NegBinEmissionScheme(0.5, 100.0).sampler()
        val sMed = NegBinEmissionScheme(10.0, 30.0).sampler()
        val sHigh = NegBinEmissionScheme(40.0, 1.0).sampler()

        val emissions = IntArray(10000) { sLow.asInt + 1 } +
                IntArray(1000) { sMed.asInt + 1 } +
                IntArray(100) { sHigh.asInt + 1 }
        emissions.shuffle()
        val (means0, failures0, _, _) = NegBinUtil.guessByData(
            emissions, 3
        )
        assertEquals(1.7, means0[0], .1)
        assertEquals(6.2, means0[1], .1)
        assertEquals(23.6, means0[2], .1)
        assertEquals(0.1, failures0[0], .1)
        assertEquals(0.8, failures0[1], .1)
        assertEquals(19.8, failures0[2], .1)
    }

}