package org.jetbrains.bio.span.statistics.hmm

import org.junit.Assert
import org.junit.Test

class AbstractNBZHMMTest {
    @Test
    fun testSNR() {
        val snrs = DoubleArray(5) { FreeNBZHMM.multiStartSignalToNoise(it) }
        Assert.assertTrue(doubleArrayOf(20.0, 40.0, 10.0, 80.0, 5.0).contentEquals(snrs))
    }
}