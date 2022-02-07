package org.jetbrains.bio.span.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.span.fit.ZLH
import org.jetbrains.bio.span.fit.ZLHID
import org.jetbrains.bio.span.statistics.emission.NegBinUtil.guessByData
import org.jetbrains.bio.span.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.hmm.MLConstrainedHMM
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.stochastic
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory

/**
 * A HMM with multidimensional Negative Binomial emissions and
 * dedicated singular zero emission, which allows specifying
 * equivalence classes between state-dimension pairs.
 *
 * The model comes with two fitters [Fitter] specially tailored for
 * detecting the presence and the difference of ChIP-Seq enrichment.
 *
 * @author Alexey Dievsky
 * @author Sergei Lebedev
 * @author Oleg Shpynov
 * @since 20/05/14
 */
class ConstrainedNBZHMM(
    stateDimensionEmissionMap: Array<IntArray>,
    means: DoubleArray,
    failures: DoubleArray
) : MLConstrainedHMM(
    stateDimensionEmissionMap,
    F64Array.stochastic(stateDimensionEmissionMap.size),
    F64Array.stochastic(stateDimensionEmissionMap.size, stateDimensionEmissionMap.size)
) {

    private val zeroEmissionScheme = ConstantIntegerEmissionScheme(0)
    private val negBinEmissionSchemes = Array(numEmissionSchemes - 1) {
        NegBinEmissionScheme(means[it], failures[it])
    }

    val means: F64Array
        get() = F64Array(numEmissionSchemes - 1) {
            negBinEmissionSchemes[it].mean
        }

    val failures: F64Array
        get() = F64Array(numEmissionSchemes - 1) {
            negBinEmissionSchemes[it].failures
        }

    val successProbabilities: F64Array
        get() = F64Array(numEmissionSchemes - 1) {
            negBinEmissionSchemes[it].successProbability
        }

    operator fun get(e: Int) = getEmissionScheme(e)

    operator fun set(e: Int, scheme: NegBinEmissionScheme) {
        if (e == 0) throw IllegalArgumentException()
        negBinEmissionSchemes[e - 1] = scheme
    }

    override fun getEmissionScheme(e: Int): IntegerEmissionScheme {
        return if (e == 0) zeroEmissionScheme else negBinEmissionSchemes[e - 1]
    }

    override fun toString(): String = toStringHelper()
        .add("means", means)
        .add("failures", failures)
        .toString()

    companion object {
        @Suppress("unused", "MayBeConstant")
        @Transient
        @JvmField
        val VERSION = 2

        private val LOG = LoggerFactory.getLogger(ConstrainedNBZHMM::class.java)

        /**
         * The fitter for detecting the presence of ChIP-Seq enrichment based on one or several
         * ChIP-Seq replicates. It assumes [ZLH] state set.
         */
        fun fitter(numReplicates: Int): Fitter<ConstrainedNBZHMM> {
            return object : Fitter<ConstrainedNBZHMM> {
                override fun guess(
                    preprocessed: Preprocessed<DataFrame>,
                    title: String,
                    threshold: Double,
                    maxIter: Int
                ): ConstrainedNBZHMM {
                    val df = preprocessed.get()
                    val means = DoubleArray(numReplicates * 2)
                    val failures = DoubleArray(numReplicates * 2)
                    for (d in 0 until numReplicates) {
                        // Filter out 0s, since they are covered by dedicated ZERO state
                        val data = df.sliceAsInt(df.labels[d]).filter { it != 0 }.sorted()
                        check(data.isNotEmpty()) {
                            "Model can't be trained on empty coverage, exiting."
                        }
                        LOG.debug("Replicate $d")
                        val (meansD, failuresD) = guessByData(data, 2)
                        means[d] = meansD[0]
                        means[d + numReplicates] = meansD[1]
                        failures[d] = failuresD[0]
                        failures[d + numReplicates] = failuresD[1]
                    }

                    return ConstrainedNBZHMM(ZLH.constraintMap(numReplicates), means, failures)
                }

                override fun fit(
                    preprocessed: List<Preprocessed<DataFrame>>,
                    title: String,
                    threshold: Double,
                    maxIter: Int,
                ) = super.fit(preprocessed, title, threshold, maxIter).apply {
                    flipStatesIfNecessary(numReplicates)
                }
            }
        }

        /**
         * The fitter for detecting differential ChIP-Seq enrichment based on one or several
         * ChIP-Seq replicates for each track.
         * It assumes [ZLHID] state set.
         */
        fun fitter(numReplicates1: Int, numReplicates2: Int): Fitter<ConstrainedNBZHMM> {
            return object : Fitter<ConstrainedNBZHMM> {
                override fun guess(
                    preprocessed: Preprocessed<DataFrame>,
                    title: String,
                    threshold: Double,
                    maxIter: Int
                ) = guess(listOf(preprocessed), title, threshold, maxIter)

                override fun guess(
                    preprocessed: List<Preprocessed<DataFrame>>,
                    title: String,
                    threshold: Double,
                    maxIter: Int
                ): ConstrainedNBZHMM {
                    val df = DataFrame.rowBind(preprocessed.map { it.get() }.toTypedArray())
                    val means = DoubleArray((numReplicates1 + numReplicates2) * 2)
                    val failures = DoubleArray((numReplicates1 + numReplicates2) * 2)
                    for (d1 in 0 until numReplicates1) {
                        // Filter out 0s, since they are covered by dedicated ZERO state
                        val data1 = df.sliceAsInt(df.labels[d1]).filter { it != 0 }.sorted()
                        check(data1.isNotEmpty()) {
                            "Model can't be trained on empty coverage (track $d1), exiting."
                        }
                        LOG.debug("Replicate $d1")
                        val (meansD1, failuresD1) = guessByData(data1, 2)
                        means[d1] = meansD1[0]
                        means[d1 + numReplicates1] = meansD1[1]
                        failures[d1] = failuresD1[0]
                        failures[d1 + numReplicates1] = failuresD1[1]
                    }
                    for (d2 in 0 until numReplicates2) {
                        // Filter out 0s, since they are covered by dedicated ZERO state
                        val data2 = df.sliceAsInt(df.labels[d2 + numReplicates1]).filter { it != 0 }.sorted()
                        check(data2.isNotEmpty()) {
                            "Model can't be trained on empty coverage (track ${d2 + numReplicates1}), exiting."
                        }
                        LOG.debug("Replicate $d2")
                        val (meansD2, failuresD2) = guessByData(data2, 2)
                        means[d2 + numReplicates1 * 2] = meansD2[0]
                        means[d2 + numReplicates2 + numReplicates1 * 2] = meansD2[1]
                        failures[d2 + numReplicates1 * 2] = failuresD2[0]
                        failures[d2 + numReplicates2 + numReplicates1 * 2] = failuresD2[1]
                    }
                    return ConstrainedNBZHMM(ZLHID.constraintMap(numReplicates1, numReplicates2), means, failures)
                }

                override fun fit(
                    preprocessed: List<Preprocessed<DataFrame>>,
                    title: String,
                    threshold: Double,
                    maxIter: Int
                ) = super.fit(preprocessed, title, threshold, maxIter).apply {
                    flipStatesIfNecessary(numReplicates1, numReplicates2)
                }
            }
        }

        /**
         * Rearranges the transition and prior probabilities so that the provided states flip
         */
        internal fun ConstrainedNBZHMM.probabilityFlip(state1: Int, state2: Int, stateNum: Int) {
            for (s in 0 until stateNum) {
                val tmp = this.logTransitionProbabilities[s, state1]
                this.logTransitionProbabilities[s, state1] = this.logTransitionProbabilities[s, state2]
                this.logTransitionProbabilities[s, state2] = tmp
            }
            for (s in 0 until stateNum) {
                val tmp = this.logTransitionProbabilities[state1, s]
                this.logTransitionProbabilities[state1, s] = this.logTransitionProbabilities[state2, s]
                this.logTransitionProbabilities[state2, s] = tmp
            }
            val tmp = this.logPriorProbabilities[state1]
            this.logPriorProbabilities[state1] = this.logPriorProbabilities[state2]
            this.logPriorProbabilities[state2] = tmp
        }

        /**
         * Flip states in case when states with HIGH get lower mean than LOW
         */
        internal fun ConstrainedNBZHMM.flipStatesIfNecessary(tracks: Int) {
            val means = this.means
            val ps = this.successProbabilities
            val switchNeeded = (0 until tracks).filter { means[it] > means[it + tracks] && ps[it] > ps[it + tracks] }
            val switchNotNeeded = (0 until tracks).filter { means[it] < means[it + tracks] || ps[it] < ps[it + tracks] }
            if (switchNeeded.isNotEmpty() && switchNotNeeded.isNotEmpty()) {
                LOG.error("Irrecoverable fitting error")
                LOG.error("means: $means")
                LOG.error("ps: $ps")
                LOG.error(
                    "track(s) " + switchNeeded.joinToString(transform = Int::toString)
                            + " contradict track(s) " + switchNotNeeded.joinToString(transform = Int::toString)
                )
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
        internal fun ConstrainedNBZHMM.flipStatesIfNecessary(tracks1: Int, tracks2: Int) {
            val means = means
            val ps = successProbabilities
            val switchNeeded1 =
                (0 until tracks1).filter { means[it] > means[it + tracks1] && ps[it] > ps[it + tracks1] }
            val switchNotNeeded1 =
                (0 until tracks1).filter { means[it] < means[it + tracks1] && ps[it] < ps[it + tracks1] }
            val switchNeeded2 = (2 * tracks1 until 2 * tracks1 + tracks2)
                .filter { means[it] > means[it + tracks2] && ps[it] > ps[it + tracks2] }
            val switchNotNeeded2 = (2 * tracks1 until 2 * tracks1 + tracks2)
                .filter { means[it] < means[it + tracks2] && ps[it] < ps[it + tracks2] }
            if (switchNeeded1.isNotEmpty() && switchNotNeeded1.isNotEmpty()) {
                LOG.error("Irrecoverable fitting error")
                LOG.error("means: $means")
                LOG.error("ps: $ps")
                LOG.error(
                    "track(s) " + switchNeeded1.joinToString(transform = Int::toString)
                            + " contradict track(s) " + switchNotNeeded1.joinToString(transform = Int::toString)
                )
                throw IllegalStateException("Irrecoverable fitting error")
            }
            if (switchNeeded2.isNotEmpty() && switchNotNeeded2.isNotEmpty()) {
                LOG.error("Irrecoverable fitting error")
                LOG.error("means: $means")
                LOG.error("ps: $ps")
                LOG.error(
                    "track(s) " + switchNeeded2.joinToString(transform = Int::toString)
                            + " contradict track(s) " + switchNotNeeded2.joinToString(transform = Int::toString)
                )
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

    }
}
