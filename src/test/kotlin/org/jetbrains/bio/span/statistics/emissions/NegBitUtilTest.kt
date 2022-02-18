package org.jetbrains.bio.span.statistics.emissions

import org.jetbrains.bio.span.statistics.emission.NegBinUtil
import org.junit.Assert
import org.junit.Test
import kotlin.test.assertEquals

class NegBitUtilTest {
    @Test
    fun testSNR() {
        val snrs = DoubleArray(5) { NegBinUtil.multiStartSignalToNoise(it) }
        Assert.assertTrue(doubleArrayOf(20.0, 40.0, 10.0, 80.0, 5.0).contentEquals(snrs))
    }

    @Test
    fun testGuessByData() {
        val emissions = IntArray(20) {1} + IntArray(20) {10} + IntArray(20) {200}
        val (means0, failures0) = NegBinUtil.guessByData(
            emissions, 2, 0
        )
        assertEquals(3.51, means0[0], 1e-2)
        assertEquals(70.33, means0[1], 1e-2)
        assertEquals(11.64, failures0[0], 1e-2)
        assertEquals(77.36, failures0[1], 1e-2)

        val (means3, failures3) = NegBinUtil.guessByData(
            emissions, 2, 3
        )
        assertEquals(0.87, means3[0], 1e-2)
        assertEquals(70.33, means3[1], 1e-2)
        assertEquals(46.59, failures3[0], 1e-2)
        assertEquals(77.36, failures3[1], 1e-2)

        val (means4, failures4) = NegBinUtil.guessByData(
            emissions, 2, 4
        )
        assertEquals(14.06, means4[0], 1e-2)
        assertEquals(70.33, means4[1], 1e-2)
        assertEquals(15.47, failures4[0], 1e-2)
        assertEquals(77.36, failures4[1], 1e-2)

    }
}