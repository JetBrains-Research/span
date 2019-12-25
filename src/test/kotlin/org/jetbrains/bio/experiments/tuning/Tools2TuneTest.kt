package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.experiments.tuning.Tool2Tune.Companion.combineParams
import org.junit.Test
import kotlin.test.assertEquals

class Tools2TuneTest {

    @Test
    fun testCombineParamsEmpty() {
        assertEquals(emptyList(), combineParams())
    }

    @Test
    fun testCombineParams() {
        val combinations = combineParams(listOf(1, 2, 3, 4, 5) to 2).map { it.toList() }
        assertEquals("[[3], [4], [2], [5], [1]]", combinations.toString())
    }

    @Test
    fun testCombineParamsDouble() {
        val combinations = combineParams(listOf(-1, 0, 1) to 1, listOf(10, 20, 30) to 0).map { it.toList() }
        assertEquals("[[0, 10], [0, 20], [0, 30], [1, 10], [1, 20], [1, 30], [-1, 10], [-1, 20], [-1, 30]]",
                combinations.toString())
    }

}