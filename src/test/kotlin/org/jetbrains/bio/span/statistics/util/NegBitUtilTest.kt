package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.distribution.Sampling
import org.junit.Test
import kotlin.test.assertEquals

class NegBitUtilTest {

    @Test
    fun testGuessByData() {
        Sampling.RANDOM_DATA_GENERATOR.reSeed(100)
        val s1 = NegBinEmissionScheme(40.0, 1.0).sampler()
        val s2 = NegBinEmissionScheme(0.5, 100.0).sampler()

        val emissions = IntArray(100) { s1.asInt + 1 } + IntArray(10000) { s2.asInt + 1 }
        val (means0, failures0, _, _) = NegBinUtil.guessByData(
            emissions, 2
        )
        assertEquals(1.3, means0[0], .1)
        assertEquals(9.3, means0[1], .1)
        assertEquals(0.1, failures0[0], .1)
        assertEquals(5.1, failures0[1], .1)
    }
}