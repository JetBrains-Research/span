package org.jetbrains.bio.span.fit

import org.jetbrains.bio.span.SpanCLA
import org.jetbrains.bio.span.peaks.MultipleTesting
import org.jetbrains.bio.viktor.F64Array
import java.lang.reflect.Modifier
import kotlin.math.ln

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

    const val SPAN_AUTOCORRELATION_MAX_SHIFT = 50

    const val SPAN_AUTOCORRELATION_CHECKPOINT = 10

    const val SPAN_FRAGMENTATION_MAX_GAP = 50

    const val SPAN_FRAGMENTATION_GAP_CHECKPOINT = 20

    // Minimal coefficient between variance and mean of Negative Binomials
    const val SPAN_NB_VAR_MEAN_MULTIPLIER = 1.1

    // Fraction top scores used for HMM high state initialization,
    // used to estimate signal-to-noise
    // Useful for narrow marks even with high noise,
    // however too small value leads to fragmentation of broad marks
    const val SPAN_INITIAL_SCORES_HIGH = 0.1

    // Fraction low scores used for HMM low state initialization
    // used to estimate signal-to-noise
    const val SPAN_INITIAL_SCORES_LOW = 0.6

    // Keep SPAN noise state at least estimated noise * this multiplier,
    // generally low mean ~1 and low value prevents fragmentation in broad marks
    const val SPAN_NOISE_MULTIPLIER = 0.4

    // k27me3 tweaked initialization priors from low signal-to-noise ratio ChIP-seq
    val SPAN_NB2ZHMM_PRIORS = F64Array.of(0.75, 0.24, 0.01)

    val SPAN_NB2ZHMM_TRANSITIONS = listOf(
        doubleArrayOf(0.75, 0.249, 0.001),
        doubleArrayOf(0.2, 0.78, 0.02),
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
     * Sensitivity background configures threshold of candidate enriched bins, allowing for bigger number of candidates
     * to be checked by enrichment vs. control track. The less value the more candidates are presented for check.
     * Potential drawback - mixing two adjacent peaks into a single one, or too much noise.
     */
    val SPAN_DEFAULT_SENSITIVITY = ln(1e-3)

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

     // Minimal distance between candidates
    const val SPAN_DEFAULT_GAP = 3

    /**
     * Clipping allows to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    const val SPAN_SIGNAL_CLIP = 0.5

    const val SPAN_LENGTH_CLIP = 0.5

    /**
     * Array of steps used for reducing the range by [SPAN_CLIP_STEPS] from both sides while increasing score.
     */
    val SPAN_CLIP_STEPS = intArrayOf(10, 20, 50, 100, 200, 500, 1000)

    fun printSpanConstants() {
        for (field in SpanConstants::class.java.declaredFields) {
            // Ignore singleton and fields, configured with params
            if (field.name in setOf(
                    "INSTANCE",
                    "SPAN_DEFAULT_BIN", "SPAN_DEFAULT_FDR",
                    "SPAN_DEFAULT_SENSITIVITY", "SPAN_DEFAULT_GAP",
                    "SPAN_FIT_THRESHOLD", "SPAN_FIT_MAX_ITERATIONS"
                )) {
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