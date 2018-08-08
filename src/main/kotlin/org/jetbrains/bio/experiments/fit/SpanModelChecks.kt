package org.jetbrains.bio.experiments.fit

import org.apache.log4j.Logger
import org.jetbrains.bio.experiments.tuning.SPAN
import org.jetbrains.bio.statistics.NegBinEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM


private val LOG = Logger.getLogger(SPAN::class.java)

internal fun MLFreeNBHMM.probabilityFlip(state1: Int, state2: Int) {
    for (i in 0..2) {
        val tmp = this.logTransitionProbabilities[i, state1]
        this.logTransitionProbabilities[i, state1] = this.logTransitionProbabilities[i, state2]
        this.logTransitionProbabilities[i, state2] = tmp
    }
    for (j in 0..2) {
        val tmp = this.logTransitionProbabilities[state1, j]
        this.logTransitionProbabilities[state1, j] = this.logTransitionProbabilities[state2, j]
        this.logTransitionProbabilities[state2, j] = tmp
    }
    val tmp = this.logPriorProbabilities[state1]
    this.logPriorProbabilities[state1] = this.logPriorProbabilities[state2]
    this.logPriorProbabilities[state2] = tmp
}


/**
 * Rearranges the transition and prior probabilities so that the provided states flip
 */
internal fun MLConstrainedNBHMM.probabilityFlip(state1: Int, state2: Int, stateNum: Int) {
    for (i in 0 until stateNum) {
        val tmp = this.logTransitionProbabilities[i, state1]
        this.logTransitionProbabilities[i, state1] = this.logTransitionProbabilities[i, state2]
        this.logTransitionProbabilities[i, state2] = tmp
    }
    for (j in 0 until stateNum) {
        val tmp = this.logTransitionProbabilities[state1, j]
        this.logTransitionProbabilities[state1, j] = this.logTransitionProbabilities[state2, j]
        this.logTransitionProbabilities[state2, j] = tmp
    }
    val tmp = this.logPriorProbabilities[state1]
    this.logPriorProbabilities[state1] = this.logPriorProbabilities[state2]
    this.logPriorProbabilities[state2] = tmp
}

/**
 * Flip states in case when states with HIGH get lower mean than LOW
 */
internal fun MLFreeNBHMM.flipStatesIfNecessary() {
    val lowScheme = this[1] as NegBinEmissionScheme
    val highScheme = this[2] as NegBinEmissionScheme
    val meanLow = means[0]
    val meanHigh = means[1]
    val pLow = lowScheme.successProbability
    val pHigh = highScheme.successProbability
    val meanFlipped = meanLow > meanHigh
    val pFlipped = pLow > pHigh
    if (meanFlipped) {
        LOG.warn("After fitting the model, mean emission in LOW state ($meanLow) is higher than " +
                "mean emission in HIGH state ($meanHigh).")
    }
    if (pFlipped) {
        LOG.warn("After fitting the model, emission's parameter p in LOW state ($pLow) is higher than " +
                "emission's parameter p in HIGH state ($pHigh).")
    }
    if (meanFlipped && pFlipped) {
        LOG.warn("This usually indicates that the states were flipped during fitting. We will now flip them back.")
        this[2] = lowScheme
        this[1] = highScheme
        this.probabilityFlip(1, 2)
    } else if (meanFlipped || pFlipped) {
        LOG.warn("This is generally harmless, but could indicate low quality of data.")
    }
}

/**
 * Flip states in case when states with HIGH get lower mean than LOW
 */
internal fun MLConstrainedNBHMM.flipStatesIfNecessary(tracks: Int) {
    val means = this.means
    val ps = this.successProbabilities
    val switchNeeded = (0 until tracks).filter { means[it] > means[it + tracks] && ps[it] > ps[it + tracks] }
    val switchNotNeeded = (0 until tracks).filter { means[it] < means[it + tracks] || ps[it] < ps[it + tracks] }
    if (switchNeeded.isNotEmpty() && switchNotNeeded.isNotEmpty()) {
        LOG.error("Irrecoverable fitting error")
        LOG.error("means: " + means.toString())
        LOG.error("ps: " + ps.toString())
        LOG.error("track(s) " + switchNeeded.joinToString(transform = Int::toString)
                + " contradict track(s) " + switchNotNeeded.joinToString(transform = Int::toString))
        throw IllegalStateException("Irrecoverable fitting error")
    }
    if (switchNeeded.isNotEmpty()) {
        // flip LOW(1) <-> HIGH(2)
        for (i in 0 until tracks) {
            val tmp = this[i + 1]
            this[i + 1] = this[i + tracks + 1] as NegBinEmissionScheme
            this[i + tracks + 1] = tmp as NegBinEmissionScheme
        }
        this.probabilityFlip(1, 2, 3)
    }
    this.updateTransients()
}

/**
 * Flip states in case when states with HIGH get lower mean than LOW
 */
internal fun MLConstrainedNBHMM.flipStatesIfNecessary(tracks1: Int, tracks2: Int) {
    val means = means
    val ps = successProbabilities
    val switchNeeded1 = (0 until tracks1).filter { means[it] > means[it + tracks1] && ps[it] > ps[it + tracks1] }
    val switchNotNeeded1 = (0 until tracks1).filter { means[it] < means[it + tracks1] && ps[it] < ps[it + tracks1] }
    val switchNeeded2 = (2 * tracks1 until 2 * tracks1 + tracks2)
            .filter { means[it] > means[it + tracks2] && ps[it] > ps[it + tracks2] }
    val switchNotNeeded2 = (2 * tracks1 until 2 * tracks1 + tracks2)
            .filter { means[it] < means[it + tracks2] && ps[it] < ps[it + tracks2] }
    if (switchNeeded1.isNotEmpty() && switchNotNeeded1.isNotEmpty()) {
        LOG.error("Irrecoverable fitting error")
        LOG.error("means: " + means.toString())
        LOG.error("ps: " + ps.toString())
        LOG.error("track(s) " + switchNeeded1.joinToString(transform = Int::toString)
                + " contradict track(s) " + switchNotNeeded1.joinToString(transform = Int::toString))
        throw IllegalStateException("Irrecoverable fitting error")
    }
    if (switchNeeded2.isNotEmpty() && switchNotNeeded2.isNotEmpty()) {
        LOG.error("Irrecoverable fitting error")
        LOG.error("means: " + means.toString())
        LOG.error("ps: " + ps.toString())
        LOG.error("track(s) " + switchNeeded2.joinToString(transform = Int::toString)
                + " contradict track(s) " + switchNotNeeded2.joinToString(transform = Int::toString))
        throw IllegalStateException("Irrecoverable fitting error")
    }
    if (switchNeeded1.isNotEmpty()) {
        // flip LOW(0) <-> DECREASING(2) and HIGH(3) <-> INCREASING(1)
        for (i in 0 until tracks1) {
            val tmp = this[i + 1]
            this[i + 1] = this[i + tracks1 + 1] as NegBinEmissionScheme
            this[i + tracks1 + 1] = tmp as NegBinEmissionScheme
        }
        this.probabilityFlip(0, 2, 5)
        this.probabilityFlip(1, 3, 5)
    }
    if (switchNeeded2.isNotEmpty()) {
        // flip LOW(0) <-> INCREASING(1) and HIGH(3) <-> DECREASING(2)
        for (i in 2 * tracks1 until 2 * tracks1 + tracks2) {
            val tmp = this[i + 1]
            this[i + 1] = this[i + tracks2 + 1] as NegBinEmissionScheme
            this[i + tracks2 + 1] = tmp as NegBinEmissionScheme
        }
        this.probabilityFlip(0, 1, 5)
        this.probabilityFlip(2, 3, 5)
    }
    updateTransients()
}