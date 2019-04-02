package org.jetbrains.bio.statistics.hmm

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.*
import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.emission.ConstantIntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.emission.NegBinEmissionScheme
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.statistics.state.ZLHID
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
class MLConstrainedNBHMM(stateDimensionEmissionMap: Array<IntArray>,
                         means: DoubleArray, failures: DoubleArray)
    : MLConstrainedHMM(stateDimensionEmissionMap,
        F64Array.stochastic(stateDimensionEmissionMap.size),
        F64Array.stochastic(stateDimensionEmissionMap.size,
                stateDimensionEmissionMap.size)) {

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
        @Transient
        @JvmField
        val VERSION = 1

        fun preprocessor() = Preprocessors.identity()

        /**
         * The fitter for detecting the presence of ChIP-Seq enrichment based on one or several
         * ChIP-Seq replicates. It assumes [ZLH] state set.
         */
        fun fitter(numReplicates: Int): Fitter<MLConstrainedNBHMM> {
            return object : Fitter<MLConstrainedNBHMM> {
                override fun guess(preprocessed: Preprocessed<DataFrame>, title: String,
                                   threshold: Double, maxIter: Int): MLConstrainedNBHMM {
                    val df = preprocessed.get()
                    val meanCoverage = DoubleArray(numReplicates)
                    val means = DoubleArray(numReplicates * 2)
                    val failures = DoubleArray(numReplicates * 2)
                    for (d in 0 until numReplicates) {
                        val values = df.sliceAsInt(df.labels[d])
                        meanCoverage[d] = values.average()
                        // Extreme initial states configuration
                        means[d] = meanCoverage[d] / Math.sqrt(MLFreeNBHMM.SIGNAL_TO_NOISE_RATIO)
                        means[d + numReplicates] = meanCoverage[d] * Math.sqrt(MLFreeNBHMM.SIGNAL_TO_NOISE_RATIO)
                        val sd = values.standardDeviation()
                        failures[d] = NegativeBinomialDistribution
                                .estimateFailuresUsingMoments(meanCoverage[d], sd * sd)
                        failures[d + numReplicates] = failures[d]
                    }

                    return MLConstrainedNBHMM(ZLH.constraintMap(numReplicates),
                            means, failures)
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
                                   threshold: Double, maxIter: Int): MLConstrainedNBHMM =
                        guess(listOf(preprocessed), title, threshold, maxIter)

                override fun guess(preprocessed: List<Preprocessed<DataFrame>>, title: String,
                                   threshold: Double, maxIter: Int): MLConstrainedNBHMM {
                    val df = DataFrame.rowBind(preprocessed.map { it.get() }.toTypedArray())
                    val means = DoubleArray((numReplicates1 + numReplicates2) * 2)
                    val failures = DoubleArray((numReplicates1 + numReplicates2) * 2)
                    for (d1 in 0 until numReplicates1) {
                        val values = df.sliceAsInt(df.labels[d1])
                        val meanCoverage = values.average()
                        check(meanCoverage != 0.0) {
                            "Model can't be trained on empty coverage " +
                                    "(track $d1), exiting."
                        }
                        // Extreme initial states configuration
                        means[d1] = meanCoverage / Math.sqrt(MLFreeNBHMM.SIGNAL_TO_NOISE_RATIO)
                        means[d1 + numReplicates1] = meanCoverage * Math.sqrt(MLFreeNBHMM.SIGNAL_TO_NOISE_RATIO)
                        val sd = values.standardDeviation()
                        failures[d1] = NegativeBinomialDistribution
                                .estimateFailuresUsingMoments(meanCoverage, sd * sd)
                        failures[d1 + numReplicates1] = failures[d1]
                    }
                    for (d2 in 0 until numReplicates2) {
                        val values = df.sliceAsInt(df.labels[d2 + numReplicates1])
                        val meanCoverage = values.average()
                        check(meanCoverage != 0.0) {
                            "Model can't be trained on empty coverage " +
                                    "(track ${d2 + numReplicates1}), exiting."
                        }
                        means[d2 + numReplicates1 * 2] =
                                meanCoverage / Math.sqrt(MLFreeNBHMM.SIGNAL_TO_NOISE_RATIO)
                        means[d2 + numReplicates2 + numReplicates1 * 2] =
                                meanCoverage * Math.sqrt(MLFreeNBHMM.SIGNAL_TO_NOISE_RATIO)
                        val sd = values.standardDeviation()
                        failures[d2 + numReplicates1 * 2] = NegativeBinomialDistribution
                                .estimateFailuresUsingMoments(meanCoverage, sd * sd)
                        failures[d2 + numReplicates2 + numReplicates1 * 2] = failures[d2]
                    }
                    return MLConstrainedNBHMM(
                            ZLHID.constraintMap(numReplicates1, numReplicates2),
                            means, failures)
                }
            }
        }
    }
}
