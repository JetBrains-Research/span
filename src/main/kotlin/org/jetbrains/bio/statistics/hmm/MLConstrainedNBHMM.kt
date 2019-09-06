package org.jetbrains.bio.statistics.hmm

import org.apache.log4j.Logger
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.flipStatesIfNecessary
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.standardDeviation
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.statistics.state.ZLHID
import org.jetbrains.bio.statistics.stochastic
import org.jetbrains.bio.viktor.F64Array

/**
 * A HMM with multidimensional Negative Binomial emissions and
 * dedicated singular zero emission, which allows to specify
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
class MLConstrainedNBHMM(
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
        val VERSION = 1

        private val LOG = Logger.getLogger(MLConstrainedNBHMM::class.java)

        /**
         * The fitter for detecting the presence of ChIP-Seq enrichment based on one or several
         * ChIP-Seq replicates. It assumes [ZLH] state set.
         */
        fun fitter(numReplicates: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> {
                override fun guess(preprocessed: Preprocessed<DataFrame>, title: String,
                        threshold: Double, maxIter: Int, attempt: Int): MLConstrainedNBHMM {
                    val df = preprocessed.get()
                    val meanCoverage = DoubleArray(numReplicates)
                    val means = DoubleArray(numReplicates * 2)
                    val failures = DoubleArray(numReplicates * 2)
                    for (d in 0 until numReplicates) {
                        // Filter out 0s, since they are covered by dedicated ZERO state
                        val values = df.sliceAsInt(df.labels[d]).filter { it != 0 }.toIntArray()
                        check(values.isNotEmpty()) {
                            "Model can't be trained on empty coverage, exiting."
                        }
                        val mean = values.average()
                        val sd = values.standardDeviation()
                        val snr = MLFreeNBHMM.signalToNoise(attempt)
                        val meanLow = mean / Math.sqrt(snr)
                        val meanHigh = mean * Math.sqrt(snr)
                        val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, sd * sd)
                        LOG.debug("Guess $attempt emissions mean $mean\tsd $sd")
                        LOG.debug("Guess $attempt init meanLow $meanLow\tmeanHigh $meanHigh\tfailures $fs")
                        meanCoverage[d] = mean
                        means[d] = meanLow
                        means[d + numReplicates] = meanHigh
                        failures[d] = fs
                        failures[d + numReplicates] = fs
                    }

                    return MLConstrainedNBHMM(ZLH.constraintMap(numReplicates), means, failures)
                }

                override fun fit(
                        preprocessed: List<Preprocessed<DataFrame>>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): MLConstrainedNBHMM {
                    return super.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                        flipStatesIfNecessary(numReplicates)
                    }
                }
            }
        }

        /**
         * The fitter for detecting differential ChIP-Seq enrichment based on one or several
         * ChIP-Seq replicates for each track.
         * It assumes [ZLHID] state set.
         */
        fun fitter(numReplicates1: Int, numReplicates2: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> {
                override fun guess(preprocessed: Preprocessed<DataFrame>, title: String,
                        threshold: Double, maxIter: Int, attempt: Int): MLConstrainedNBHMM =
                        guess(listOf(preprocessed), title, threshold, maxIter, attempt)

                override fun guess(preprocessed: List<Preprocessed<DataFrame>>, title: String,
                        threshold: Double, maxIter: Int, attempt: Int): MLConstrainedNBHMM {
                    val df = DataFrame.rowBind(preprocessed.map { it.get() }.toTypedArray())
                    val means = DoubleArray((numReplicates1 + numReplicates2) * 2)
                    val failures = DoubleArray((numReplicates1 + numReplicates2) * 2)
                    for (d1 in 0 until numReplicates1) {
                        // Filter out 0s, since they are covered by dedicated ZERO state
                        val values = df.sliceAsInt(df.labels[d1]).filter { it != 0 }.toIntArray()
                        check(values.isNotEmpty()) {
                            "Model can't be trained on empty coverage " +
                                    "(track $d1), exiting."
                        }
                        val mean = values.average()
                        val sd = values.standardDeviation()
                        val snr = MLFreeNBHMM.signalToNoise(attempt)
                        val meanLow = mean / Math.sqrt(snr)
                        val meanHigh = mean * Math.sqrt(snr)
                        val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, sd * sd)
                        LOG.debug("Guess1 $attempt emissions mean $mean\tsd $sd")
                        LOG.debug("Guess1 $attempt init meanLow $meanLow\tmeanHigh $meanHigh\tfailures $fs")
                        means[d1] = meanLow
                        means[d1 + numReplicates1] = meanHigh
                        failures[d1] = fs
                        failures[d1 + numReplicates1] = fs
                    }
                    for (d2 in 0 until numReplicates2) {
                        // Filter out 0s, since they are covered by dedicated ZERO state
                        val values = df.sliceAsInt(df.labels[d2 + numReplicates1]).filter { it != 0 }.toIntArray()
                        check(values.isNotEmpty()) {
                            "Model can't be trained on empty coverage " +
                                    "(track ${d2 + numReplicates1}), exiting."
                        }
                        val mean = values.average()
                        val sd = values.standardDeviation()
                        val snr = MLFreeNBHMM.signalToNoise(attempt)
                        val meanLow = mean / Math.sqrt(snr)
                        val meanHigh = mean * Math.sqrt(snr)
                        val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, sd * sd)
                        LOG.debug("Guess2 $attempt emissions mean $mean\tsd $sd")
                        LOG.debug("Guess2 $attempt init meanLow $meanLow\tmeanHigh $meanHigh\tfailures $fs")
                        means[d2 + numReplicates1 * 2] = meanLow
                        means[d2 + numReplicates2 + numReplicates1 * 2] = meanHigh
                        failures[d2 + numReplicates1 * 2] = fs
                        failures[d2 + numReplicates2 + numReplicates1 * 2] = fs
                    }
                    return MLConstrainedNBHMM(
                        ZLHID.constraintMap(numReplicates1, numReplicates2),
                        means, failures)
                }

                override fun fit(
                        preprocessed: List<Preprocessed<DataFrame>>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): MLConstrainedNBHMM {
                    return super.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                        flipStatesIfNecessary(numReplicates1, numReplicates2)
                    }
                }
            }
        }
    }
}
