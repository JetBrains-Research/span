package org.jetbrains.bio.span.statistics.emissions

import org.jetbrains.bio.span.statistics.emission.NegBinUtil
import org.junit.Assert
import org.junit.Test
import kotlin.test.assertEquals

class NegBitUtilTest {
    @Test
    fun testSNR() {
        val snrs = DoubleArray(5) { NegBinUtil.multiStartSignalToNoise(it) }
        Assert.assertTrue(doubleArrayOf(10.0, 20.0, 5.0, 40.0, 2.5).contentEquals(snrs))
    }

    @Test
    fun testGuessByData() {
        val emissions = IntArray(20) {1} + IntArray(20) {10} + IntArray(20) {200}
        val (means0, failures0) = NegBinUtil.guessByData(
            emissions, 2, 0
        )
        assertEquals(70.33, means0[0], 1e-2)
        assertEquals(703.33, means0[1], 1e-2)
        assertEquals(77.36, failures0[0], 1e-2)
        assertEquals(773.66, failures0[1], 1e-2)
    }
}