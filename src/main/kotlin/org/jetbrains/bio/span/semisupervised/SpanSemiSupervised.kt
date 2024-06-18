package org.jetbrains.bio.span.semisupervised

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FDR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.semisupervised.LocationLabel.Companion.computeErrors
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.MultitaskProgress
import org.jetbrains.bio.util.await
import java.util.concurrent.Callable

object SpanSemiSupervised {

    private val SPAN_FDRS = listOf(
        0.1, SPAN_DEFAULT_FDR, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10
    )

    val SPAN_BACKGROUND_SENSITIVITY_VARIANTS =
        doubleArrayOf(10.0, 5.0, 2.0, 1.2, 1.0, 0.8, 0.5, 0.2, 0.1, 0.05, 0.001, 1e-4, 1e-6, 1e-8).sorted()

    val SPAN_GAPS_VARIANTS =
        doubleArrayOf(0.0, SPAN_DEFAULT_GAP, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0).sorted()

    val PARAMETERS =
        SPAN_FDRS.sorted().flatMap { fdr ->
            SPAN_BACKGROUND_SENSITIVITY_VARIANTS.sorted().flatMap { sensitivity ->
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
        parameters: List<Triple<Double, Double, Double>>,
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

        val tasks = parameters.mapIndexed { index, (fdr, sensitivity, clip) ->
            Callable {
                cancellableState.checkCanceled()
                val peaksOnLabeledGenomeQuery =
                    ModelToPeaks.getPeaks(
                        results, labeledGenomeQuery, fdr, sensitivity, clip,
                        CancellableState.current()
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