package org.jetbrains.bio.statistics.state

import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.junit.Assert.assertArrayEquals
import org.junit.Test
import kotlin.test.assertEquals
import kotlin.test.assertFalse
import kotlin.test.assertTrue

/**
 * Created by Aleksei Dievskii on 04.07.2018.
 */


class ChipSeqTest {

    @Test
    fun testConstraints() {
        val zhlConstraints = ZLH.constraintMap(3)
        assertAllEqual(zhlConstraints[ZLH.Z.ordinal],
                "Emission schemes for Z state should be equal")
        assertNotAllEqual(zhlConstraints[ZLH.L.ordinal],
                "Emission schemes for L state should be different")
        assertNotAllEqual(zhlConstraints[ZLH.H.ordinal],
                "Emission schemes for H state should be different")
        val zlhidConstraints = ZLHID.constraintMap(3, 2)
        assertAllEqual(zlhidConstraints[ZLHID.Z.ordinal],
                "Emission schemes for Z state should be equal")
        assertNotAllEqual(zlhidConstraints[ZLHID.L.ordinal],
                "Emission schemes for L state should be different")
        assertNotAllEqual(zlhidConstraints[ZLHID.I.ordinal],
                "Emission schemes for I state should be different")
        assertNotAllEqual(zlhidConstraints[ZLHID.D.ordinal],
                "Emission schemes for D state should be different")
        assertNotAllEqual(zlhidConstraints[ZLHID.H.ordinal],
                "Emission schemes for H state should be different")
    }

    @Test
    fun testModel() {
        val zhlConstraints = ZLH.constraintMap(1)
        val free = MLFreeNBHMM(0.5, 5.0, 1.0)
        assertIs(free[zhlConstraints[ZLH.Z.ordinal][0]], ConstantIntegerEmissionScheme::class.java)
        assertIs(free[zhlConstraints[ZLH.L.ordinal][0]], NegBinEmissionScheme::class.java)
        assertIs(free[zhlConstraints[ZLH.H.ordinal][0]], NegBinEmissionScheme::class.java)

        val zlhidConstraints = ZLHID.constraintMap(3, 2)
        val constrained = MLConstrainedNBHMM(zlhidConstraints,
                                             DoubleArray(10) { it.toDouble() },
                                             DoubleArray(10) { it.toDouble() })
        zlhidConstraints[ZLHID.Z.ordinal]
                .forEach { assertIs(constrained[it], ConstantIntegerEmissionScheme::class.java) }
        zlhidConstraints[ZLHID.L.ordinal]
                .forEach { assertIs(constrained[it], NegBinEmissionScheme::class.java) }
        zlhidConstraints[ZLHID.I.ordinal]
                .forEach { assertIs(constrained[it], NegBinEmissionScheme::class.java) }
        zlhidConstraints[ZLHID.D.ordinal]
                .forEach { assertIs(constrained[it], NegBinEmissionScheme::class.java) }
        zlhidConstraints[ZLHID.H.ordinal]
                .forEach { assertIs(constrained[it], NegBinEmissionScheme::class.java) }
    }

    @Test
    fun testZLHIDSemantics() {
        val zlhidConstraints = ZLHID.constraintMap(3, 2)
        assertArrayEquals("L and I emissions differ for first group",
                          zlhidConstraints[ZLHID.L.ordinal].sliceArray(0 until 3),
                          zlhidConstraints[ZLHID.I.ordinal].sliceArray(0 until 3))
        assertArrayEquals("D and H emissions differ for first group",
                          zlhidConstraints[ZLHID.D.ordinal].sliceArray(0 until 3),
                          zlhidConstraints[ZLHID.H.ordinal].sliceArray(0 until 3))
        assertArrayEquals("L and D emissions differ for second group",
                          zlhidConstraints[ZLHID.L.ordinal].sliceArray(3 until 5),
                          zlhidConstraints[ZLHID.D.ordinal].sliceArray(3 until 5))
        assertArrayEquals("I and H emissions differ for second group",
                          zlhidConstraints[ZLHID.I.ordinal].sliceArray(3 until 5),
                          zlhidConstraints[ZLHID.H.ordinal].sliceArray(3 until 5))
    }

    companion object {
        fun assertAllEqual(array: IntArray, message: String) {
            if (array.isEmpty()) return
            val item = array.first()
            array.forEach { assertEquals(item, it, message) }
        }

        fun assertNotAllEqual(array: IntArray, message: String) {
            assertFalse(array.isEmpty(), "$message (array is empty)")
            val item = array.first()
            assertFalse(array.all { item == it }, "$message (all elements are equal to $item)")
        }

        fun assertIs(actual: Any, expected: Class<out Any>) {
            assertTrue(expected.isInstance(actual),
                       "Expected ${expected.simpleName}, got ${actual::class.java.simpleName}.")
        }
    }
}