package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.withTempFile
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
import java.util.stream.Collectors
import kotlin.test.assertEquals

/**
 * @author Oleg Shpynov
 * @since 04/04/2018.
 */

class CoverageFitInformationTest {
    @get:Rule
    var expectedEx = ExpectedException.none()

    val gq = GenomeQuery("to1")
    val chr2 = gq.get()[1]

    @Test
    fun checkBinSize() {
        expectedEx.expect(IllegalStateException::class.java)
        expectedEx.expectMessage("Wrong bin size, expected: 100, got: 50")
        CoverageFitInformation("foo", 100, "hg19", LinkedHashMap()).checkBinSize(50)
    }

    @Test
    fun checkBuild() {
        expectedEx.expect(IllegalStateException::class.java)
        expectedEx.expectMessage("Wrong genome build, expected: hg19, got: to2")
        CoverageFitInformation("foo", 100, "hg19", LinkedHashMap()).checkGenomeQuery(GenomeQuery("to2"))
    }

    @Test
    fun checkGenomeQuery() {
        expectedEx.expect(IllegalStateException::class.java)
        expectedEx.expectMessage("Wrong chromosomes, expected: [chr1, chr2, chr3, chrX], got: [chr1]")
        CoverageFitInformation.of("foo", 100, gq).checkGenomeQuery(GenomeQuery("to1", "chr1"))
    }

    @Test
    fun checkGenomeQueryOrder() {
        CoverageFitInformation.of("foo", 100, GenomeQuery("to1", "chr1", "chr2"))
                .checkGenomeQuery(GenomeQuery("to1", "chr2", "chr1"))
    }



    val JSON = """{
  "description": "foobar",
  "bin_size": 200,
  "build": "to1",
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrX": 1000000
  }
}
"""

    @Test
    fun checkOf() {
        val of = CoverageFitInformation.of("foobar", 200, gq)
        assertEquals(listOf("chr1", "chr2", "chr3", "chrX"), of.chromosomesSizes.keys.toList())
    }

    @Test
    fun checkSave() {
        val info = CoverageFitInformation.of("foobar", 200, gq)
        withTempFile("foo", ".tar") { path ->
            info.save(path)
            val output = path.bufferedReader().lines().collect(Collectors.joining("\n"))
            assertEquals(JSON.trim(), output.trim())
        }
        assertEquals(listOf("chr1", "chr2", "chr3", "chrX"), info.chromosomesSizes.keys.toList())
    }

    @Test
    fun checkLoad() {
        val info = CoverageFitInformation.of("foobar", 200, gq)
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use { it.write(JSON) }
            assertEquals(info, CoverageFitInformation.load(path))
        }
    }

    @Test
    fun checkIndices() {
        val info = CoverageFitInformation.of("foobar", 200, gq)
        assertEquals(50000 to 55000, info.getChromosomesIndices(gq, chr2))
    }

    @Test
    fun checkOffsets() {
        val info = CoverageFitInformation.of("foobar", 200, gq)
        assertEquals(listOf(0, 200, 400, 600, 800), info.offsets(chr2).take(5))
        assertEquals(5000, chr2.length / 200)
        assertEquals(5000, info.offsets(chr2).size)
    }
}

