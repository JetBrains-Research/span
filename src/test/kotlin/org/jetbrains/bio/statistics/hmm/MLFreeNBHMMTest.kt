package org.jetbrains.bio.statistics.hmm

import org.junit.Assert
import org.junit.Test

class MLFreeNBHMMTest {
    @Test
    fun testSNR() {
        val snrs = DoubleArray(5) { MLFreeNBHMM.multiStartSignalToNoise(it) }
        Assert.assertTrue(doubleArrayOf(20.0, 40.0, 10.0, 80.0, 5.0).contentEquals(snrs))
    }
}