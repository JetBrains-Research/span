package org.jetbrains.bio.span.fit

import org.jetbrains.bio.viktor.F64Array

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
     * The default step size used for beta values in the SPAN calculation.
     * Beta value is used to minimize correlation between treatment and control track.
     * Step is used to iterate in the [0, 1] interval.
     */
    const val SPAN_DEFAULT_BETA_STEP = 0.01

    /**
     * Default signal-to-noise ratio used for HMM model states initialization
     */
    const val SPAN_DEFAULT_SIGNAL_TO_NOISE = 20.0

    val SPAN_DEFAULT_NB2ZHMM_PRIORS = F64Array.of(0.6, 0.399, 0.001)

    val SPAN_DEFAULT_NB2ZHMM_TRANSITIONS_ = listOf(
        doubleArrayOf(0.5, 0.499, 0.001),
        doubleArrayOf(0.399, 0.6, 0.001),
        doubleArrayOf(0.05, 0.15, 0.8))

    /**
     * SPAN_FIT_THRESHOLD is a constant variable that represents the threshold value used
     * in the SPAN peak-fitting algorithm, relative delta between the model log likelihood
     */
    const val SPAN_FIT_THRESHOLD = 1e-3

    /**
     * The maximum number of iterations to perform during the fitting process in the SPAN algorithm.
     */
    const val SPAN_FIT_MAX_ITERATIONS = 10

    /**
     * Sensitivity background configures threshold of candidate enriched bins, allowing for bigger number of candidates
     * to be checked by enrichment vs. control track. The less value the more candidates are presented for check.
     * Potential drawback - mixing two adjacent peaks into a single one, or too much noise.
     */
    const val SPAN_DEFAULT_BACKGROUND_SENSITIVITY = 0.2

    /**
     * Fraction of top significant bins to be used for peak score estimation.
     * It should be robust wrt appending blocks of low significance,
     * so take into account top N% bins into block, otherwise we'll get fdr-blinking peaks, i.e.,
     * peaks which are present for stronger fdr, but missing for more relaxed settings
     */
    const val SPAN_SCORE_BLOCKS = 0.5

    /**
     * Minimal distance between score blocks
     */
    const val SPAN_SCORE_BLOCKS_DISTANCE = 0.05

    /**
     * Clipping allows to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    const val SPAN_DEFAULT_CLIP = 0.4

    /**
     * Array of steps used for reducing the range by [SPAN_CLIP_STEPS] from both sides while increasing score.
     */
    val SPAN_CLIP_STEPS = intArrayOf(10, 20, 50, 100, 200, 500, 1000)

    /**
     * Threshold to ignore chromosome with different coverage than genome-wide.
     * Useful for contigs, where coverage may be much higher than generally across the genome.
     */
    const val SPAN_MAXIMUM_COVERAGE_DELTA = 0.01
}