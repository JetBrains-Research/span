package org.jetbrains.bio.span.fit

import org.jetbrains.bio.span.peaks.ModelToPeaks.relaxedLogFdr

object SpanConstants {
    /**
     * Default bin size used for binned centered reads coverage aggregation.
     */
    const val SPAN_DEFAULT_BIN = 200

    /**
     * The SPAN_DEFAULT_FDR variable represents the default value for the FDR parameter.
     * The FDR is a measure of the expected proportion of false positive peaks in the peak calling results.
     */
    const val SPAN_DEFAULT_FDR = 0.05

    /**
     * Default gap value used to merge peaks, particularly useful for broad histone marks peak calling.
     */
    const val SPAN_DEFAULT_GAP = 3

    /**
     * The default step size used for beta values in the SPAN calculation.
     * Beta value is used to minimize correlation between treatment and control track.
     * Step is used to iterate in the [0, 1] interval.
     */
    const val SPAN_DEFAULT_BETA_STEP = 0.01

    /**
     * SPAN_FIT_THRESHOLD is a constant variable that represents the threshold value used
     * in the SPAN peak fitting algorithm.
     */
    const val SPAN_FIT_THRESHOLD = 1.0

    /**
     * The maximum number of iterations to perform during the fitting process in the SPAN algorithm.
     */
    const val SPAN_FIT_MAX_ITERATIONS = 20

    /**
     * Peak calling is done in two runs - with relaxed settings and strict user-provided strict FDR.
     * This value is used to compute relaxed settings, see [relaxedLogFdr] for details.
     */
    const val SPAN_RELAX_POWER_DEFAULT = 0.5

    /**
     * Clipping allows to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    const val SPAN_DEFAULT_CLIP = true

    /**
     * Array of steps used for reducing the range by [SPAN_CLIP_STEPS] from both sides while increasing score.
     * Used together with [SPAN_MAX_CLIPPED_FRACTION].
     */
    val SPAN_CLIP_STEPS = intArrayOf(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000)

    /**
     * The maximum fraction of length for clipped peaks, from left and from right.
     */
    const val SPAN_MAX_CLIPPED_FRACTION = 0.25

    /**
     * Maximum threshold value used for clipping peaks by coverage.
     * Defines max score of clipped region as average_noise + MAX_CLIPPED_THRESHOLD * (average_signal - average_noise)
     */
    const val SPAN_MAX_CLIPPED_THRESHOLD = 0.75

}