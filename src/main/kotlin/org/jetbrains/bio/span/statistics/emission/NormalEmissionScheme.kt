package org.jetbrains.bio.span.statistics.emission

import org.jetbrains.bio.span.standardDeviation
import org.jetbrains.bio.statistics.distribution.NormalIntDistribution
import org.jetbrains.bio.statistics.emission.IntegerEmissionScheme
import org.jetbrains.bio.statistics.gson.NotDirectlyDeserializable
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.LoggerFactory
import java.util.function.IntSupplier

class NormalEmissionScheme(mean: Double, variance: Double) : IntegerEmissionScheme, NotDirectlyDeserializable {

    override val degreesOfFreedom: Int = 2

    var mean: Double = mean
        internal set

    var variance: Double = variance
        internal set

    init {
        require(variance >= 0) { "Variance must be non-negative, got $variance" }
        updateTransients()
    }

    @Transient
    private var backNormalIntDistribution = NormalIntDistribution(mean, variance)

    override fun updateTransients() {
        backNormalIntDistribution = NormalIntDistribution(mean, variance)
    }


    override fun sampler(): IntSupplier {
        return IntSupplier { backNormalIntDistribution.sample() }
    }

    override fun logProbability(value: Int): Double = when {
        mean.isNaN() || variance.isNaN() || value == Integer.MIN_VALUE || value == Integer.MAX_VALUE -> Double.NEGATIVE_INFINITY
        else -> backNormalIntDistribution.logProbability(value)
    }

    override fun update(sample: IntArray, weights: F64Array) {
        mean = weights.dot(sample) / weights.sum()
        val sd = DoubleArray(sample.size) { sample[it] * weights[it] }.standardDeviation()
        variance = sd * sd
        updateTransients()
        LOG.debug("Mean $mean\tVariance $variance")
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(NormalEmissionScheme::class.java)
    }
}

