package org.jetbrains.bio.span.fit

import org.jetbrains.bio.statistics.model.DifferenceType


/**
 * Zero inflated enrichment status.
 * Used in SPAN enrichment analysis.
 */
enum class ZLH {
    Z,  // ZERO
    L,  // LOW
    H;  // HIGH

    companion object {
        // XXX unlike other implementations this one produces transposed map.
        @JvmStatic
        fun constraintMap(numReplicates: Int): Array<IntArray> {
            val res = Array(3) { IntArray(numReplicates) }
            for (d in 0 until numReplicates) {
                res[1][d] = d + 1
                res[2][d] = d + numReplicates + 1
            }
            return res
        }
    }
}

/**
 * Zero inflated INCREASED / DECREASED state.
 * Used in SPAN differential analysis.
 */
enum class ZLHID : DifferenceType<ZLH> {
    L,  // LOW
    I,  // INCREASED
    D,  // DECREASED
    H,  // HIGH
    Z;  // ZERO

    override val isDifferent: Boolean get() = this in different()

    override fun first() = when (this) {
        Z -> ZLH.Z
        L, I -> ZLH.L
        D, H -> ZLH.H
    }

    override fun second() = when (this) {
        Z -> ZLH.Z
        L, D -> ZLH.L
        I, H -> ZLH.H
    }

    companion object {

        fun constraintMap(numReplicates1: Int, numReplicates2: Int): Array<IntArray> {
            val res = Array(5) { IntArray(numReplicates1 + numReplicates2) }
            for (d in 0 until numReplicates1) {
                res[0][d] = d + 1
                res[1][d] = d + 1
                res[2][d] = d + 1 + numReplicates1
                res[3][d] = d + 1 + numReplicates1
            }
            for (d in 0 until numReplicates2) {
                res[0][d + numReplicates1] = d + 1 + 2 * numReplicates1
                res[1][d + numReplicates1] = d + 1 + numReplicates2 + 2 * numReplicates1
                res[2][d + numReplicates1] = d + 1 + 2 * numReplicates1
                res[3][d + numReplicates1] = d + 1 + numReplicates2 + 2 * numReplicates1
            }
            return res
        }

        fun different() = setOf(I, D)

        @JvmStatic
        fun same() = setOf(Z, L, H)
    }
}
