package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.util.glob
import org.jetbrains.bio.util.withTempDirectory
import org.jetbrains.bio.util.withTempFile
import org.junit.Test
import kotlin.test.assertEquals

/**
 * @author Oleg Shpynov
 * @date 2019-03-07
 */
class DataConfigExtTest {
    @Test
    fun testFillData() {
        withConfig { config ->
            withTempDirectory("symlinks") { dir ->
                val list = config.fillData(dir, null)
                val expected = listOf(
                        "test0_input.bam",
                        "test0_rep1_H3K27ac.bam",
                        "test0_rep2_H3K27ac.bam",
                        "test1_input.bam",
                        "test1_r0_H3K9me3.bam",
                        "test1_rep1_H3K4me3.bam",
                        "test2_input.bam",
                        "test2_r0_H3K4me3.bam")
                assertEquals(expected,
                        dir.glob("*.bam").map { it.fileName.toString() }.sorted())
                assertEquals(expected, list.map { it.fileName.toString() }.sorted())
            }
        }
    }

    @Test
    fun testFillDataMark() {
        withConfig { config ->
            withTempDirectory("symlinks") { dir ->
                config.fillData(dir, "H3K4me3", addInput = false)
                assertEquals(listOf(
                        "test1_rep1_H3K4me3.bam",
                        "test2_r0_H3K4me3.bam"),
                        dir.glob("*.bam").map { it.fileName.toString() }.sorted())
            }
        }
    }

    private fun withConfig(block: (DataConfig) -> Unit) {
        withTempFile("test1", ".bam") { p11 ->
            withTempFile("test1", ".bam") { p12 ->
                withTempFile("test2", ".bam") { p2 ->
                    withTempFile("input", ".bam") { input ->
                        val config = DataConfig.load("""genome: to1
tracks:
   Input:
      test0:
      - $input
      test1:
      - $input
      test2:
      - $input
   H3K27ac:
      test0:
        rep1:
          path: $p11
          meta: true
        rep2:
          path: $p12
          meta: true
   H3K4me3:
      test1:
        rep1: $p11
        rep2:
           path: $p12
           failed: true
      test2:
      - $p2
   H3K9me3:
      test1:
      - $p2
""".reader(), "test")
                        block(config)
                    }
                }
            }
        }
    }

}