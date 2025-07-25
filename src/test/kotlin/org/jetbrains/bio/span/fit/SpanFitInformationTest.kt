package org.jetbrains.bio.span.fit

import com.google.gson.JsonParseException
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.withTempFile
import org.junit.Assert.assertEquals
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
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
        SpanAnalyzeFitInformation(
            "hg19", emptyList(), null, emptyList(), FixedFragment(100), true, 200, LinkedHashMap()
        ).checkGenome(Genome["to1"])
    }

    @Test
    fun checkGenomeQueryOrder() {
        SpanAnalyzeFitInformation(
            GenomeQuery(Genome["to1"], "chr1", "chr2"), emptyList(), null, emptyList(), FixedFragment(100), true, 200
        ).checkGenome(Genome["to1"])
    }


    @Test
    fun checkOf() {
        val of = SpanAnalyzeFitInformation(
            gq, emptyList(), null, emptyList(), FixedFragment(100), true, 200
        )
        assertEquals(listOf("chr1", "chr2", "chr3", "chrM", "chrX"), of.chromosomesSizes.keys.toList())
    }

    @Test
    fun checkSave() {
        withTempFile("treatment", ".bam") { t ->
            withTempFile("control", ".bam") { c ->
                val info = SpanAnalyzeFitInformation(
                    gq, listOf(SpanDataPaths(t, c)), null, listOf("treatment_control"),
                    FixedFragment(100), false, 200
                )
                withTempFile("foo", ".tar") { path ->
                    info.save(path)
                    // Escape Windows separators here
                    assertEquals(
                        """{
  "build": "to1",
  "paths": [
    {
      "treatment": "${t.toString().replace("\\", "\\\\")}",
      "control": "${c.toString().replace("\\", "\\\\")}"
    }
  ],
  "labels": [
    "treatment_control"
  ],
  "fragment": "100",
  "unique": false,
  "bin_size": 200,
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrM": 10000,
    "chrX": 1000000
  },
  "fit.information.fqn": "org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation",
  "version": 5
}""".trim().lines(), path.bufferedReader().lines().collect(Collectors.toList())
                    )
                }
                assertEquals(listOf("chr1", "chr2", "chr3", "chrM", "chrX"), info.chromosomesSizes.keys.toList())
            }
        }
    }


    @Test
    fun checkWrongVersion() {
        expectedEx.expect(JsonParseException::class.java)
        expectedEx.expectMessage("expects '${SpanAnalyzeFitInformation.VERSION}' version, but got '100500'")
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use {
                it.write(
                    """{
  "build": "to1",
  "paths": [],
  "labels": [],
  "fragment": "auto",
  "bin_size": 200,
  "ie6_compatibility": false,
  "enable_quantum_optimization": true,
  "chromosomes_sizes": {},
  "fit.information.fqn": "org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation",  
  "version": 100500
}"""
                )
            }
            SpanFitInformation.load<SpanFitInformation>(path)
        }
    }

    @Test
    fun checkWrongFqn() {
        withTempFile("foo", ".tar") { path ->
            expectedEx.expect(JsonParseException::class.java)
            expectedEx.expectMessage(
                "Cannot load class org.jetbrains.bio.span.fit.Span100500FitInformation"
            )
            path.bufferedWriter().use {
                it.write(
                    """{
  "build": "to1",
  "paths": [
    {
      "treatment": "path_to_file",
      "control": "path_to_control"
    }
  ],
  "mapability_path": "mapability.bigWig",
  "labels": [
    "treatment_control"
  ],
  "fragment": 42,
  "unique": false,
  "bin_size": 200,
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrM": 1000000,
    "chrX": 1000000
  },
  "fit.information.fqn": "org.jetbrains.bio.span.fit.Span100500FitInformation",
  "version": 3
}"""
                )
            }
            SpanFitInformation.load<SpanFitInformation>(path)
        }
    }

    @Test
    fun checkNoVersion() {
        withTempFile("foo", ".tar") { path ->
            expectedEx.expect(JsonParseException::class.java)
            expectedEx.expectMessage("Version field (version) is missing")
            path.bufferedWriter().use {
                it.write(
                    """{
  "foo": "bar",
  "baz": [],
  "fit.information.fqn": "org.jetbrains.bio.span.fit.SpanFitInformation"
}"""
                )
            }
            SpanFitInformation.load<SpanFitInformation>(path)
        }
    }

    @Test
    fun checkNoFqn() {
        withTempFile("foo", ".tar") { path ->
            expectedEx.expect(JsonParseException::class.java)
            expectedEx.expectMessage("Class name (fit.information.fqn) is missing")
            path.bufferedWriter().use {
                it.write(
                    """{
    "build": "to1",
  "paths": [
    {
      "treatment": "path_to_file",
      "control": "path_to_control"
    }
  ],
  "mapability_path": "mapability.bigWig",
  "labels": [
    "treatment_control"
  ],
  "fragment": 42,
  "unique": false,
  "bin_size": 200,
  "chromosomes_sizes": {
    "chr1": 10000000,
    "chr2": 1000000,
    "chr3": 1000000,
    "chrM": 1000000,
    "chrX": 1000000
  },
  "version": 5
}"""
                )
            }
            SpanFitInformation.load<SpanFitInformation>(path)
        }
    }

    @Test
    fun checkIndices() {
        val info = SpanAnalyzeFitInformation(
            gq, emptyList(), null, emptyList(), FixedFragment(100), true, 200
        )
        assertEquals(50000 to 55000, info.getChromosomesIndices(chr2))
    }

    @Test
    fun checkOffsets() {
        val info = SpanAnalyzeFitInformation(
            gq, emptyList(), null, emptyList(), FixedFragment(100), true, 200
        )
        assertEquals(listOf(0, 200, 400, 600, 800), info.offsets(chr2).take(5))
        assertEquals(5000, chr2.length / 200)
        assertEquals(5000, info.offsets(chr2).size)
    }
}