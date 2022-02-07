package org.jetbrains.bio.span.statistics.emission

import org.jetbrains.bio.statistics.distribution.NegativeBinomialDistribution
import org.jetbrains.bio.statistics.standardDeviation
import org.slf4j.LoggerFactory
import kotlin.math.max

object NegBinUtil {

    private val LOG = LoggerFactory.getLogger(NegBinUtil::class.java)

    /**
     * Guess initial approximation for mixture / HMM given data for several Negative Binomial distributions.
     */
    fun guessByData(
        sortedData: List<Int>,
        n: Int
    ): Pair<DoubleArray, DoubleArray> {
        val means = DoubleArray(n)
        val failures = DoubleArray(n)
        (0 until n).forEach { i ->
            val subList = sortedData.subList(sortedData.size / n * i, sortedData.size / n * (i + 1)).toIntArray()
            val mean = max(0.1 * (i + 1), subList.average())
            val sd = subList.standardDeviation()
            // NegativeBinomialDistribution requires variance greater than mean, tweak variance if required.
            // Otherwise, failures will be set to +Inf and won't be updated.
            val fs = NegativeBinomialDistribution.estimateFailuresUsingMoments(mean, max(1.1 * mean, sd * sd))
            means[i] = mean
            failures[i] = fs
        }
        LOG.debug("Guess emissions means ${means.joinToString(",")}, failures ${failures.joinToString(",")}")
        return means to failures
    }
}