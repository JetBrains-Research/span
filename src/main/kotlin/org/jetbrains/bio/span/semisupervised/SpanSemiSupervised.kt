package org.jetbrains.bio.span.semisupervised

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.span.fit.SpanConstants
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FDR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_HARD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_SPEED
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_LIGHT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MIN_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.SpanModelToPeaks
import org.jetbrains.bio.span.semisupervised.LocationLabel.Companion.computeErrors
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.MultitaskProgress
import org.jetbrains.bio.util.await
import java.util.concurrent.Callable

object SpanSemiSupervised {

    private val SPAN_FDRS = listOf(
        0.1, SPAN_DEFAULT_FDR, 0.01, 1e-3, 1e-4, 1e-6, 1e-8, 1e-10, 1e-20
    )

    val SPAN_SENSITIVITY_LOG_VARIANTS: DoubleArray =
        SpanModelToPeaks.linSpace(-100.0, SPAN_MIN_SENSITIVITY, 50)

    val SPAN_GAPS_VARIANTS =
        intArrayOf(0, 2, 5, 10, 20, 50).sorted()

    val PARAMETERS =
        SPAN_FDRS.sorted().flatMap { fdr ->
            SPAN_SENSITIVITY_LOG_VARIANTS.sorted().flatMap { sensitivity ->
                SPAN_GAPS_VARIANTS.map { gap ->
                    Triple(fdr, sensitivity, gap)
                }
            }
        }

    fun tuneParameters(
        results: SpanFitResults,
        genomeQuery: GenomeQuery,
        labels: List<LocationLabel>,
        id: String,
        parameters: List<Triple<Double, Double, Int>>,
        cancellableState: CancellableState
    ): Pair<List<LabelErrors>, Int> {
        val labeledGenomeQuery = GenomeQuery(
            genomeQuery.genome,
            *labels.map { it.location.chromosome.name }.distinct().toTypedArray()
        )
        // Parallelism is OK here:
        // 1. getPeaks creates BitterSet per each parameters' combination of size
        // ~ 3*10^9 / 200bp / 8 / 1024 / 1024 ~ 2MB for human genome
        // 2. List.parallelStream()....collect(Collectors.toList()) guarantees the same order as in original list.
        // Order is important!
        val labelErrorsGrid = Array<LabelErrors?>(parameters.size) { null }

        val tasks = parameters.mapIndexed { index, (fdr, sensitivity, gap) ->
            Callable {
                cancellableState.checkCanceled()
                val peaksOnLabeledGenomeQuery =
                    SpanModelToPeaks.getPeaks(
                        results, labeledGenomeQuery, fdr,
                        SpanConstants.SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION,
                        sensitivity, gap, false,
                        SPAN_DEFAULT_FRAGMENTATION_LIGHT,
                        SPAN_DEFAULT_FRAGMENTATION_HARD,
                        SPAN_DEFAULT_FRAGMENTATION_SPEED,
                        SPAN_DEFAULT_CLIP_MAX_SIGNAL,
                        cancellableState = CancellableState.current()
                    )
                labelErrorsGrid[index] = computeErrors(
                    labels,
                    LocationsMergingList.create(
                        labeledGenomeQuery,
                        peaksOnLabeledGenomeQuery.toList().map { it.location }.iterator()
                    )
                )
                MultitaskProgress.reportTask(id)
            }
        }
        tasks.await(parallel = true)
        val minTotalError = labelErrorsGrid.minOf { it!!.error() }
        // Parameters should return desired order for each tool
        return labelErrorsGrid.map { it!! } to
                parameters.indices.first { labelErrorsGrid[it]!!.error() == minTotalError }
    }

}