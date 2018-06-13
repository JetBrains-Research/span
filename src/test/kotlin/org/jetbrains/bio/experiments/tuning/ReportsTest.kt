package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.withTempDirectory
import org.junit.Test
import kotlin.test.assertEquals

/**
 * @author Oleg Shpynov
 * @since 29/03/2018.
 */
class ReportsTest {

    @Test
    fun testDonor() {
        withTempDirectory("foo") { d ->
            assertEquals("OD14", donor(d / "OD_OD14_K27me3_hg19.bam"))
        }
    }

    @Test
    fun testParams() {
        withTempDirectory("foo") { d ->
            assertEquals("200_1.0E-4_5", params(d / "YD_YD1_H3K27ac_200_1.0E-4_5_peaks.bed", "H3K27ac"))
            assertEquals("broad1.0E-6", params(d / "OD_OD14_H3K27ac_broad1.0E-6_peaks.broadPeak", "H3K27ac"))
            assertEquals("W200-G600-FDR1.0E-8", params(d / "OD_OD8_H3K27me3-W200-G600-FDR1.0E-8-island.bed", "H3K27me3"))
        }
    }

}

