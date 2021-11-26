package org.jetbrains.bio.span

import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals

/**
 * @author Oleg Shpynov
 * @date 04/07/2017
 */

class PeaksTest {

    val genomeQuery = GenomeQuery(Genome["to1"])
    val chromosome = genomeQuery.get().first()
    val chromosome2 = genomeQuery.get()[1]

    @Test
    fun testCompare() {
        assertEquals(
            -1,
            Peak(chromosome, 0, 100, 5.0, 7.5, score = 7).compareTo(Peak(chromosome2, 0, 100, 5.0, 7.5, score = 7))
        )
        assertEquals(
            -1,
            Peak(chromosome, 0, 100, 5.0, 7.5, score = 7).compareTo(Peak(chromosome, 100, 200, 5.0, 7.5, score = 7))
        )
        assertEquals(
            -1,
            Peak(chromosome, 0, 50, 5.0, 7.5, score = 7).compareTo(Peak(chromosome, 0, 100, 5.0, 7.5, score = 7))
        )
        assertEquals(
            0,
            Peak(chromosome, 0, 100, 5.0, 7.5, score = 7).compareTo(Peak(chromosome, 0, 100, 5.0, 7.5, score = 7))
        )
    }


    @Test
    fun testSavePeaks() {
        val peaks = listOf(Peak(chromosome, 0, 100, 5.0, 7.5, score = 10))
        withTempFile("peaks", ".bed") { path ->
            Peak.savePeaks(peaks, path, "foo")
            assertEquals(
                "chr1\t0\t100\tfoo_1\t10\t.\t0.0\t5.0\t7.5",
                path.bufferedReader().readLines().joinToString("\n")
            )
        }
    }
}