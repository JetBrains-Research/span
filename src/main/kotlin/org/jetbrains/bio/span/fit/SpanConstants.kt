package org.jetbrains.bio.span.fit

import org.jetbrains.bio.viktor.F64Array
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

    /**
     * The default step size used for beta values in the SPAN calculation.
     * Beta value is used to minimize correlation between treatment and control track.
     * Step is used to iterate in the [0, 1] interval.
     */
    const val SPAN_DEFAULT_BETA_STEP = 0.01

    const val SPAN_AUTOCORRELATION_MAX_SHIFT = 50

    const val SPAN_AUTOCORRELATION_CHECKPOINT = 10

    const val SPAN_FRAGMENTATION_MAX_GAP = 50

    const val SPAN_FRAGMENTATION_GAP_CHECKPOINT = 20

    /**
     * Default signal-to-noise ratio used for HMM model states initialization
     * Too high values lead to not enough sensitivity for narrow mark cases
     */
    const val SPAN_DEFAULT_SIGNAL_TO_NOISE = 10.0

    /**
     * Keep SPAN state means values separated
     */
    const val SPAN_MIN_SIGNAL_TO_NOISE = 1.5

    const val SPAN_MAX_SIGNAL_TO_NOISE = 20.0

    const val SPAN_SIGNAL_TO_NOISE_PUSHBACK = 0.1

    // k27me3 tweaked initialization priors from low signal-to-noise ratio ChIP-seq
    val SPAN_DEFAULT_NB2ZHMM_PRIORS = F64Array.of(0.75, 0.24, 0.01)

    val SPAN_DEFAULT_NB2ZHMM_TRANSITIONS_ = listOf(
        doubleArrayOf(0.75, 0.249, 0.001),
        doubleArrayOf(0.2, 0.78, 0.02),
        doubleArrayOf(0.005, 0.015, 0.98))

    /**
     * SPAN_FIT_THRESHOLD is a constant variable that represents the threshold value used
     * in the SPAN peak-fitting algorithm, relative delta between the model log likelihood
     */
    const val SPAN_FIT_THRESHOLD = 1e-4

    /**
     * The maximum number of iterations to perform during the fitting process in the SPAN algorithm.
     */
    const val SPAN_FIT_MAX_ITERATIONS = 10

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
     * Fraction of top significant bins to be used for peak score estimation.
     * It should be robust wrt appending blocks of low significance,
     * so take into account top N% bins into block, otherwise we'll get fdr-blinking peaks, i.e.,
     * peaks which are present for stronger fdr, but missing for more relaxed settings
     */
    const val SPAN_SCORE_BLOCKS = 0.5

     // Minimal distance between bins to create core blocks in candidates
    const val SPAN_DEFAULT_GAP = 3
    /**
     * Clipping allows to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    const val SPAN_DEFAULT_SIGNAL_CLIP = 0.5

    const val SPAN_DEFAULT_LENGTH_CLIP = 0.5

    /**
     * Array of steps used for reducing the range by [SPAN_CLIP_STEPS] from both sides while increasing score.
     */
    val SPAN_CLIP_STEPS = intArrayOf(10, 20, 50, 100, 200, 500, 1000)

}