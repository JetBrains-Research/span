package org.jetbrains.bio.span.fit

import org.jetbrains.bio.span.SpanCLA
import org.jetbrains.bio.span.peaks.MultipleTesting
import org.jetbrains.bio.viktor.F64Array
import java.lang.reflect.Modifier
import kotlin.math.ln

/**
 * Constants used in SPAN.
 * Those, configurable from command line interface, have "DEFAULT" in their names.
 */
object SpanConstants {
    /**
     * Default bin size used for binned centered reads coverage aggregation.
     */
    const val SPAN_DEFAULT_BIN = 100

    /**
     * The SPAN_DEFAULT_FDR variable represents the default value for the FDR parameter.
     * The FDR is a measure of the expected proportion of false positive peaks in the peak calling results.
     */
    const val SPAN_DEFAULT_FDR = 0.05

    val SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION = MultipleTesting.BH

    /**
     * The default step size used for beta values in the SPAN calculation.
     * Beta value is used to minimize correlation between treatment and control track.
     * Step is used to iterate in the [0, 1] interval.
     */
    const val SPAN_BETA_STEP = 0.01

    // Max shift to compute auto correlations
    const val SPAN_AUTOCORRELATION_MAX_SHIFT = 50

    // Max gap to compute fragmentation,
    // i.e. reduction of candidate number when merging with gap
    const val SPAN_FRAGMENTATION_MAX_GAP = 50

    // Technical minimal coefficient between variance and mean of Negative Binomials
    const val SPAN_HMM_NB_VAR_MEAN_MULTIPLIER = 1.1

    // Fraction scores used for HMM signal estimation, guards for good signal-to-noise ratio
    const val SPAN_HMM_SIGNAL_ESTIMATE = 0.05

    // Minimal low state mean threshold, guards against over-peak calling and too broad peaks
    const val SPAN_HMM_LOW_THRESHOLD = 0.3

    // Technical threshold to limit mean to std, guards against artificial data without noise
    const val SPAN_HMM_MAX_MEAN_TO_STD = 5.0

    // General model priors and priors based on real data peaks footprint
    val SPAN_HMM_PRIORS = F64Array.of(0.75, 0.249, 0.001)

    val SPAN_HMM_TRANSITIONS = listOf(
        doubleArrayOf(0.75, 0.2499, 0.0001),
        doubleArrayOf(0.2, 0.798, 0.002),
        doubleArrayOf(0.005, 0.015, 0.98))

    /**
     * The threshold  used in the SPAN peak-fitting algorithm,
     * relative delta between the model log likelihood
     */
    const val SPAN_DEFAULT_FIT_THRESHOLD = 1e-4

    /**
     * The maximum number of iterations to perform during the fitting
     */
    const val SPAN_DEFAULT_FIT_MAX_ITERATIONS = 10

    /**
     * Number of points between relaxed and strict sensitivity to analyse
     */
    const val SPAN_SENSITIVITY_N = 100

    /**
     * Limit min sensitivity during candidates saturation analysis
     */
    const val SPAN_MIN_SENSITIVITY = -1e-10

    /**
     * Fraction of top significant bins to be used for peak score estimation.
     * It should be robust wrt appending blocks of low significance,
     * so take into account top N% bins into block, otherwise we'll get fdr-blinking peaks, i.e.,
     * peaks which are present for stronger fdr, but missing for more relaxed settings
     */
    const val SPAN_SCORE_BLOCKS = 0.5

    const val SPAN_SCORE_BLOCKS_GAP = 3

    const val SPAN_DEFAULT_GAP = 0

    val SPAN_DEFAULT_SENSITIVITY = ln(SPAN_DEFAULT_FDR)

    /**
     * Gap value is computed from model bins autocorrelation and fragmentation
     * These thresholds were estimated empirically from `--deep-analysis` SPAN output
     * after analysing autocorrelation_average_score vs fragmentation_average_score
     */

    const val SPAN_BROAD_AC_MIN_THRESHOLD = 0.5

    const val SPAN_BROAD_EXTRA_GAP = 10  // x default bin 100bp

    const val SPAN_FRAGMENTED_MAX_THRESHOLD = 30

    const val SPAN_FRAGMENTED_EXTRA_GAP = 20  // x default bin 100bp

    /**
     * Clipping allows to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    const val SPAN_CLIP_MAX_SIGNAL = 0.5

    const val SPAN_CLIP_MAX_LENGTH = 0.8

    val SPAN_CLIP_STEPS = doubleArrayOf(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.5, 2.0, 5.0, 10.0)

    fun printSpanConstants() {
        for (field in SpanConstants::class.java.declaredFields) {
            // Ignore singleton and fields, configured with params
            if (field.name == "INSTANCE" || "DEFAULT" in field.name) {
                continue
            }
            if (Modifier.isStatic(field.modifiers) && Modifier.isFinal(field.modifiers)) {
                field.isAccessible = true
                val value = field.get(null) // null because it's a static field
                SpanCLA.LOG.debug("${field.name}: ${pp(value)}")
            }
        }
    }

    fun pp(value: Any?): String {
        return when (value) {
            is List<*> -> "[${value.joinToString(", ") { pp(it) }}]"
            is IntArray -> value.contentToString()
            is DoubleArray -> value.contentToString()
            is Array<*> -> value.contentDeepToString()
            else -> value.toString()
        }
    }
}