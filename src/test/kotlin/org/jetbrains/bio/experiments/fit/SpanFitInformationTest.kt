package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.toPath
import org.jetbrains.bio.util.withTempFile
import org.junit.Assert.assertEquals
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
import java.nio.file.Path
import java.util.*
import java.util.stream.Collectors

/**
 * @author Oleg Shpynov
 * @since 04/04/2018.
 */

class SpanFitInformationTest {
    @get:Rule
    var expectedEx = ExpectedException.none()

    val gq = GenomeQuery(Genome["to1"])
    val chr2 = gq.get()[1]


    @Test
    fun checkWrongBuild() {
        expectedEx.expect(IllegalStateException::class.java)
        expectedEx.expectMessage("Wrong genome build, expected: hg19, got: to1")
        SpanFitInformation(
            "hg19", emptyList(), emptyList(), Optional.of(100), 200, LinkedHashMap(), 1
        ).checkGenome(Genome["to1"])
    }

    @Test
    fun checkGenomeQueryOrder() {
        SpanFitInformation(
            GenomeQuery(Genome["to1"], "chr1", "chr2"), emptyList(), emptyList(), 100, 200
        ).checkGenome(Genome["to1"])
    }


    @Test
    fun checkOf() {
        val of = SpanFitInformation(gq, emptyList(), emptyList(), 100, 200)
        assertEquals(listOf("chr1", "chr2", "chr3", "chrM", "chrX"), of.chromosomesSizes.keys.toList())
    }

    @Test
    fun checkSave() {
        withTempFile("treatment", ".bam") { t ->
            withTempFile("control", ".bam") { c ->
                val info = SpanFitInformation(gq, listOf(t to c), listOf("treatment_control"), 100, 200)
                withTempFile("foo", ".tar") { path ->
                    info.save(path)
                    // Escape Windows separators here
                    assertEquals("""{
  "build": "to1",
  "data": [
    {
      "path": "${t.toString().replace("\\", "\\\\")}",
      "control": "${c.toString().replace("\\", "\\\\")}"
    }
  ],
  "labels": [
    "treatment_control"
  ],
  "fragment": 100,
  "bin_size": 200,
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrM": 1000000,
    "chrX": 1000000
  },
  "version": 1
}""".trim().lines(), path.bufferedReader().lines().collect(Collectors.toList()))
                }
                assertEquals(listOf("chr1", "chr2", "chr3", "chrM", "chrX"), info.chromosomesSizes.keys.toList())
            }
        }
    }

    @Test
    fun checkLoad() {
        val info = SpanFitInformation(
            gq, listOf("path_to_file".toPath() to null), listOf("treatment_control"), Optional.empty(), 200
        )
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use {
                it.write("""{
  "build": "to1",
  "data": [
    {
      "path": "path_to_file"
    }
  ],
  "labels": [
    "treatment_control"
  ],
  "fragment": null,
  "bin_size": 200,
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrM": 1000000,
    "chrX": 1000000
  },
  "version": 1
}""")
            }
            assertEquals(info, SpanFitInformation.load(path))
        }
    }

    @Test
    fun checkVersion() {
        expectedEx.expect(IllegalStateException::class.java)
        expectedEx.expectMessage("Wrong version: expected: 1, got: 2")
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use {
                it.write("""{
  "build": "to1",
  "data": [],
  "labels": [],
  "fragment": null,
  "bin_size": 200,
  "chromosomes_sizes": {},
  "version": 2
}""")
            }
            SpanFitInformation.load(path)
        }
    }

    @Test
    fun checkIndices() {
        val info = SpanFitInformation(gq, emptyList(), emptyList(), 100, 200)
        assertEquals(50000 to 55000, info.getChromosomesIndices(chr2))
    }

    @Test
    fun checkOffsets() {
        val info = SpanFitInformation(gq, emptyList(), emptyList(), 100, 200)
        assertEquals(listOf(0, 200, 400, 600, 800), info.offsets(chr2).take(5))
        assertEquals(5000, chr2.length / 200)
        assertEquals(5000, info.offsets(chr2).size)
    }
}

/**
 * Simplified instance construction for tests.
 */
internal operator fun SpanFitInformation.Companion.invoke(
        genomeQuery: GenomeQuery,
        paths: List<Pair<Path, Path?>>,
        labels: List<String>,
        fragment: Int,
        binSize: Int
    ): SpanFitInformation = SpanFitInformation(genomeQuery, paths, labels, Optional.of(fragment), binSize)