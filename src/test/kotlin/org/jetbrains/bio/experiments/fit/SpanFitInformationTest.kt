package org.jetbrains.bio.experiments.fit

import com.google.gson.JsonParseException
import org.jetbrains.bio.experiments.fit.experimental.Span2FitInformation
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.util.bufferedReader
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.util.toPath
import org.jetbrains.bio.util.withTempFile
import org.junit.Assert.assertEquals
import org.junit.Rule
import org.junit.Test
import org.junit.rules.ExpectedException
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
        Span1AnalyzeFitInformation(
            "hg19", emptyList(), emptyList(), FixedFragment(100), true, 200, LinkedHashMap()
        ).checkGenome(Genome["to1"])
    }

    @Test
    fun checkGenomeQueryOrder() {
        Span1AnalyzeFitInformation(
            GenomeQuery(Genome["to1"], "chr1", "chr2"), emptyList(), emptyList(), 100, true,  200
        ).checkGenome(Genome["to1"])
    }


    @Test
    fun checkOf() {
        val of = Span1AnalyzeFitInformation(gq, emptyList(), emptyList(), 100, true, 200)
        assertEquals(listOf("chr1", "chr2", "chr3", "chrM", "chrX"), of.chromosomesSizes.keys.toList())
    }

    @Test
    fun checkSave() {
        withTempFile("treatment", ".bam") { t ->
            withTempFile("control", ".bam") { c ->
                val info = Span1AnalyzeFitInformation(
                    gq, listOf(SpanDataPaths(t, c)), listOf("treatment_control"),
                    100, false, 200
                )
                withTempFile("foo", ".tar") { path ->
                    info.save(path)
                    // Escape Windows separators here
                    assertEquals("""{
  "build": "to1",
  "data": [
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
    "chrM": 1000000,
    "chrX": 1000000
  },
  "fit.information.fqn": "org.jetbrains.bio.experiments.fit.Span1AnalyzeFitInformation",
  "version": 3
}""".trim().lines(), path.bufferedReader().lines().collect(Collectors.toList()))
                }
                assertEquals(listOf("chr1", "chr2", "chr3", "chrM", "chrX"), info.chromosomesSizes.keys.toList())
            }
        }
    }

    @Test
    fun checkLoad() {
        val info = Span2FitInformation(
                gq, SpanDataPaths("path_to_file".toPath(), "path_to_control".toPath()),
                "mapability.bigWig".toPath(), FixedFragment(42), false, 200
        )
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use {
                it.write("""{
  "build": "to1",
  "data": [
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
  "fit.information.fqn": "org.jetbrains.bio.experiments.fit.experimental.Span2FitInformation",  
  "version": 3
}""")
            }
            assertEquals(info, SpanFitInformation.load(path))
        }
    }

    @Test
    fun checkWrongVersion() {
        expectedEx.expect(JsonParseException::class.java)
        expectedEx.expectMessage("expects '${Span1AnalyzeFitInformation.VERSION}' version, but got '100500'")
        withTempFile("foo", ".tar") { path ->
            path.bufferedWriter().use {
                it.write("""{
  "build": "to1",
  "data": [],
  "labels": [],
  "fragment": "auto",
  "bin_size": 200,
  "ie6_compatibility": false,
  "enable_quantum_optimization": true,
  "chromosomes_sizes": {},
  "fit.information.fqn": "org.jetbrains.bio.experiments.fit.Span1AnalyzeFitInformation",  
  "version": 100500
}""")
            }
            SpanFitInformation.load<SpanFitInformation>(path)
        }
    }

    @Test
    fun checkWrongFqn() {
        withTempFile("foo", ".tar") { path ->
            expectedEx.expect(JsonParseException::class.java)
            expectedEx.expectMessage(
                "Cannot load class org.jetbrains.bio.experiments.fit.Span100500FitInformation"
            )
            path.bufferedWriter().use {
                it.write(
                    """{
  "build": "to1",
  "data": [
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
  "fit.information.fqn": "org.jetbrains.bio.experiments.fit.Span100500FitInformation",
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
  "fit.information.fqn": "org.jetbrains.bio.experiments.fit.SpanFitInformation"
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
  "data": [
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
  "version": 3
}"""
                )
            }
            SpanFitInformation.load<SpanFitInformation>(path)
        }
    }

    @Test
    fun checkIndices() {
        val info = Span1AnalyzeFitInformation(gq, emptyList(), emptyList(), 100, true, 200)
        assertEquals(50000 to 55000, info.getChromosomesIndices(chr2))
    }

    @Test
    fun checkOffsets() {
        val info = Span1AnalyzeFitInformation(gq, emptyList(), emptyList(), 100, true, 200)
        assertEquals(listOf(0, 200, 400, 600, 800), info.offsets(chr2).take(5))
        assertEquals(5000, chr2.length / 200)
        assertEquals(5000, info.offsets(chr2).size)
    }
}

/**
 * Simplified instance construction for tests.
 */
internal operator fun Span1AnalyzeFitInformation.Companion.invoke(
        genomeQuery: GenomeQuery,
        paths: List<SpanDataPaths>,
        labels: List<String>,
        fragment: Int,
        unique: Boolean,
        binSize: Int
): Span1AnalyzeFitInformation = Span1AnalyzeFitInformation(
    genomeQuery, paths, labels, FixedFragment(fragment), unique, binSize
)