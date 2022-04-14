package org.jetbrains.bio.span.statistics.util

import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.junit.Test
import kotlin.test.assertEquals

class NegBitUtilTest {

    @Test
    fun testGuessByData() {
        val s1 = NegBinEmissionScheme(40.0, 1.0).sampler()
        val s2 = NegBinEmissionScheme(0.5, 100.0).sampler()

        val emissions = IntArray(100) {s1.asInt + 1} + IntArray(10000) {s2.asInt + 1}
        val (means0, failures0) = NegBinUtil.guessByData(
            emissions, 2
        )
        assertEquals(1.87, means0[0], 1.0)
        assertEquals(18.66, means0[1], 1.0)
        assertEquals(2.05, failures0[0], 1.0)
        assertEquals(20.52, failures0[1], 1.0)
    }
}