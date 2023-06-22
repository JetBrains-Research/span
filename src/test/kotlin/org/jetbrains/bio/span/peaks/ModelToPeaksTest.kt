package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Range
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
import java.lang.IllegalArgumentException
import kotlin.test.assertEquals
import kotlin.test.assertTrue

/**
 * @author Oleg Shpynov
 */

class ModelToPeaksTest {

    @get:Rule
    var expectedEx = ExpectedException.none()

    @Test
    fun testEmpty() {
        val (bins, peaks, cores) = ModelToPeaks.computeBinsCoresAndPeaks(BitList(10), BitList(10), 0)
        assertEquals(0, bins.cardinality())
        assertTrue(peaks.isEmpty())
        assertTrue(cores.isEmpty())
    }

    @Test
    fun testSimple() {
        val relaxedBins = BitList(10).apply { listOf(0, 1, 2, 4, 5).forEach { set(it) } }
        val strictBins = BitList(10).apply { listOf(1, 4, 5).forEach { set(it) } }
        val (bins, peaks, cores) = ModelToPeaks.computeBinsCoresAndPeaks(relaxedBins, strictBins, 0)
        assertEquals(5, bins.cardinality())
        assertEquals(arrayListOf(Range(0, 3), Range(4, 6)), peaks)
        assertEquals(arrayListOf(listOf(Range(1, 2)), listOf(Range(4, 6))), cores)
    }


    @Test
    fun testMultipleCores() {
        val relaxedBins = BitList(10).apply { listOf(0, 1, 2, 3, 4, 5, 8, 9, 10).forEach { set(it) } }
        val strictBins = BitList(10).apply { listOf(0, 2, 3, 9).forEach { set(it) } }
        val (bins, peaks, cores) = ModelToPeaks.computeBinsCoresAndPeaks(relaxedBins, strictBins, 0)
        assertEquals(9, bins.cardinality())
        assertEquals(arrayListOf(Range(0, 6), Range(8, 10)), peaks)
        assertEquals(arrayListOf(listOf(Range(0, 1), Range(2, 4)), listOf(Range(9, 10))), cores)
    }

    @Test
    fun testGap() {
        val relaxedBins = BitList(10).apply { listOf(0, 1, 2, 3, 4, 5, 8, 9, 10).forEach { set(it) } }
        val strictBins = BitList(10).apply { listOf(0, 2, 3, 9).forEach { set(it) } }
        val (bins, peaks, cores) = ModelToPeaks.computeBinsCoresAndPeaks(relaxedBins, strictBins, 3)
        assertEquals(9, bins.cardinality())
        assertEquals(arrayListOf(Range(0, 10)), peaks)
        assertEquals(arrayListOf(listOf(Range(0, 1), Range(2, 4), Range(9, 10))), cores)
    }

    @Test
    fun emptyCores() {
        val relaxedBins = BitList(10).apply { listOf(0, 1, 2, 3, 4, 5, 8, 9, 10).forEach { set(it) } }
        val strictBins = BitList(10).apply { listOf(0, 2, 3).forEach { set(it) } }
        val (bins, peaks, cores) = ModelToPeaks.computeBinsCoresAndPeaks(relaxedBins, strictBins, 3)
        assertEquals(9, bins.cardinality())
        assertEquals(arrayListOf(Range(0, 6)), peaks)
        assertEquals(arrayListOf(listOf(Range(0, 1), Range(2, 4))), cores)
    }

    @Test
    fun differentSize() {
        expectedEx.expect(IllegalArgumentException::class.java)
        expectedEx.expectMessage("Different size")
        val relaxedBins = BitList(20).apply { listOf(0, 1, 2, 3, 4, 5, 8, 9, 10).forEach { set(it) } }
        val strictBins = BitList(10).apply { listOf(0, 2, 3).forEach { set(it) } }
        ModelToPeaks.computeBinsCoresAndPeaks(relaxedBins, strictBins, 3)
    }


    @Test
    fun strictNotRelaxes() {
        expectedEx.expect(IllegalArgumentException::class.java)
        expectedEx.expectMessage("Strict positions should be covered in relaxed")
        val relaxedBins = BitList(10).apply { listOf(1, 2, 3, 4, 5, 8, 9, 10).forEach { set(it) } }
        val strictBins = BitList(10).apply { listOf(0, 2, 3).forEach { set(it) } }
        ModelToPeaks.computeBinsCoresAndPeaks(relaxedBins, strictBins, 3)
    }
}