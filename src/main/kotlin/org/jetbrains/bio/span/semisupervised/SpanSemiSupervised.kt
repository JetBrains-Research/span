package org.jetbrains.bio.span.semisupervised

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.semisupervised.LocationLabel.Companion.computeErrors
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.MultitaskProgress
import org.jetbrains.bio.util.awaitAll
import org.jetbrains.bio.util.parallelismLevel
import java.util.concurrent.Callable
import java.util.concurrent.Executors

object SpanSemiSupervised {

    val FDRS = listOf(
        0.1,
        SpanPeakCallingExperiment.SPAN_DEFAULT_FDR, 0.01, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10
    )

    val GAPS = listOf(0, 1, 2, SpanPeakCallingExperiment.SPAN_DEFAULT_GAP, 5, 10)

    val PARAMETERS =
        FDRS.sorted().flatMap { fdr ->
            GAPS.sorted().map { gap -> fdr to gap }
        }

    private val TUNING_EXECUTOR = Executors.newWorkStealingPool(parallelismLevel())

    fun tune(
        results: SpanFitResults,
        genomeQuery: GenomeQuery,
        labels: List<LocationLabel>,
        id: String,
        parameters: List<Pair<Double, Int>>,
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
        TUNING_EXECUTOR.awaitAll(
            parameters.mapIndexed { index, (fdr, gap) ->
                Callable {
                    cancellableState.checkCanceled()
                    val peaksOnLabeledGenomeQuery =
                        ModelToPeaks.computeChromosomePeaks(
                            results, labeledGenomeQuery, fdr, gap, false,
                            CancellableState.current()
                        )
                    labelErrorsGrid[index] = computeErrors(
                        labels,
                        LocationsMergingList.create(
                            labeledGenomeQuery,
                            peaksOnLabeledGenomeQuery.map { it.location }.iterator()
                        )
                    )
                    MultitaskProgress.reportTask(id)
                }
            }

        )
        val minTotalError = labelErrorsGrid.minOf { it!!.error() }
        // Parameters should return desired order for each tool
        return labelErrorsGrid.map { it!! } to
                parameters.indices.first { labelErrorsGrid[it]!!.error() == minTotalError }
    }

}