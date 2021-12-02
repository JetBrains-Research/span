package org.jetbrains.bio.span

import org.junit.Assert
import org.junit.Test

class IslandsTest {


    // ______/\______
    //  -10  0  +10
    private val tags1 = arrayListOf<Pair<Int, Int>>().apply {
        (-10 until 0).forEach { add(it to 10 + it) }
        (0 until 10).forEach { add(it to 10 - it) }
    }


    @Test
    fun testClip1() {
        val density: (Int, Int) -> Double = { a, b ->
            tags1.filter { (x, _) -> x in a until b }.sumOf { it.second }.toDouble() / (b - a)
        }
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, 50, 0.5, density))
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, 50, 0.1, density))
        Assert.assertEquals(-9 to 10, clipIsland(-10, 10, 5, 0.25, density))
        Assert.assertEquals(-9 to 75, clipIsland(-10, 100, 25, 0.25, density))
        Assert.assertEquals(-8 to 8, clipIsland(-8, 8, 4, 0.25, density))
    }

    // ______/\/\_______
    // -20 -10 0 +10 +20
    private val tags2 = arrayListOf<Pair<Int, Int>>().apply {
        (-20 until -10).forEach { add(it to 20 + it) }
        (-10 until 0).forEach { add(it to -it) }
        (0 until 10).forEach { add(it to it) }
        (10 until 20).forEach { add(it to 20 - it) }
    }


    @Test
    fun testClip2() {
        val density: (Int, Int) -> Double = { a, b ->
            tags2.filter { (x, _) -> x in a until b }.sumOf { it.second }.toDouble() / (b - a)
        }
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, 50, 0.5, density))
        Assert.assertEquals(-16 to 17, clipIsland(-30, 30, 15, 0.5, density))
        Assert.assertEquals(-10 to 10, clipIsland(-10, 10, 5, 0.5, density))
        Assert.assertEquals(-10 to 5, clipIsland(-10, 5, 4, 0.5, density))
        Assert.assertEquals(-10 to 75, clipIsland(-10, 100, 25, 0.5, density))
        Assert.assertEquals(-8 to 8, clipIsland(-8, 8, 4, 0.5, density))
    }


    // _______/--TT--\_____
    // -20 -10 -5 5 +10 +20
    private val tags3 = arrayListOf<Pair<Int, Int>>().apply {
        (-20 until -10).forEach { add(it to 20 + it) }
        (-10 until -5).forEach { add(it to 5) }
        (-5 until 5).forEach { add(it to 20) }
        (5 until 10).forEach { add(it to 5) }
        (10 until 20).forEach { add(it to 20 - it) }
    }

    @Test
    fun testClip3() {
        val density: (Int, Int) -> Double = { a, b ->
            tags3.filter { (x, _) -> x in a until b }.sumOf { it.second }.toDouble() / (b - a)
        }
        Assert.assertEquals(-50 to 50, clipIsland(-100, 100, 50, 0.5, density))
        Assert.assertEquals(-16 to 17, clipIsland(-30, 30, 15, 0.5, density))
        Assert.assertEquals(-10 to 10, clipIsland(-10, 10, 5, 0.5, density))
        Assert.assertEquals(-10 to 5, clipIsland(-10, 5, 4, 0.5, density))
        Assert.assertEquals(-10 to 75, clipIsland(-10, 100, 25, 0.5, density))
        Assert.assertEquals(-21 to 8, clipIsland(-30, 8, 9, 0.5, density))
        Assert.assertEquals(-8 to 8, clipIsland(-8, 8, 4, 0.5, density))
    }


}