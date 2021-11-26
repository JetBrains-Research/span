package org.jetbrains.bio.span

import org.junit.Assert
import org.junit.Test

class IslandsTest {


    // ______/\______
    //  -10  0  +10
    private val scores1 = arrayListOf<Pair<Int, Int>>().apply {
        (-10 until 0).forEach { add(it to 10 + it) }
        (0 until 10).forEach { add(it to 10 - it) }
    }


    @Test
    fun testClip1() {
        val score: (Int, Int) -> Double = { a, b ->
            scores1.filter { (x, _) -> x in a until b }.sumOf { it.second }.toDouble()
        }
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, score))
        Assert.assertEquals(-10 to 10, clipIsland(-10, 10, score))
        Assert.assertEquals(-9 to 73, clipIsland(-10, 100, score))
        Assert.assertEquals(-8 to 8, clipIsland(-8, 8, score))
    }

    // ______/\/\_______
    // -20 -10 0 +10 +20
    private val scores2 = arrayListOf<Pair<Int, Int>>().apply {
        (-20 until -10).forEach { add(it to 20 + it) }
        (-10 until 0).forEach { add(it to -it) }
        (0 until 10).forEach { add(it to it) }
        (10 until 20).forEach { add(it to 20 - it) }
    }


    @Test
    fun testClip2() {
        val score: (Int, Int) -> Double = { a, b ->
            scores2.filter { (x, _) -> x in a until b }.sumOf { it.second }.toDouble()
        }
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, score))
        Assert.assertEquals(-30 to 30, clipIsland(-30, 30, score))
        Assert.assertEquals(-10 to 10, clipIsland(-10, 10, score))
        Assert.assertEquals(-10 to 5, clipIsland(-10, 5, score))
        Assert.assertEquals(-10 to 73, clipIsland(-10, 100, score))
        Assert.assertEquals(-8 to 8, clipIsland(-8, 8, score))
    }


    // _______/--TT--\_____
    // -20 -10 -5 5 +10 +20
    private val scores3 = arrayListOf<Pair<Int, Int>>().apply {
        (-20 until -10).forEach { add(it to 20 + it) }
        (-10 until -5).forEach { add(it to 5) }
        (-5 until 5).forEach { add(it to 20) }
        (5 until 10).forEach { add(it to 5) }
        (10 until 20).forEach { add(it to 20 - it) }
    }

    @Test
    fun testClip3() {
        val score: (Int, Int) -> Double = { a, b ->
            scores3.filter { (x, _) -> x in a until b }.sumOf { it.second }.toDouble()
        }
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, score))
        Assert.assertEquals(-30 to 30, clipIsland(-30, 30, score))
        Assert.assertEquals(-10 to 10, clipIsland(-10, 10, score))
        Assert.assertEquals(-10 to 5, clipIsland(-10, 5, score))
        Assert.assertEquals(-10 to 73, clipIsland(-10, 100, score))
        Assert.assertEquals(-8 to 8, clipIsland(-30, 8, score))
        Assert.assertEquals(-8 to 8, clipIsland(-8, 8, score))
    }


}