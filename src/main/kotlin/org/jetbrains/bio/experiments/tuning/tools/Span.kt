package org.jetbrains.bio.experiments.tuning.tools

import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanFitResults
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_BIN
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_FDR
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_GAP
import org.jetbrains.bio.experiments.tuning.*
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.data.Cell
import org.jetbrains.bio.genome.data.ChipSeqTarget
import org.jetbrains.bio.genome.data.DataConfig
import org.jetbrains.bio.span.peaks.*
import org.jetbrains.bio.util.*
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.Executors

object Span : Tool2Tune<Pair<Double, Int>>() {

    private val LOG = LoggerFactory.getLogger(Span::class.java)

    private val tuningExecutor = Executors.newWorkStealingPool(parallelismLevel())

    override val id = "span"
    override val suffix = "_peaks.bed"

    private val FDRS = listOf(0.1, 0.05, 1E-2, 1E-3, 1E-4, 1E-5, SPAN_DEFAULT_FDR, 1E-7, 1E-8, 1E-9, 1E-10)

    val GAPS = listOf(0, 1, 2, SPAN_DEFAULT_GAP, 10)

    override val parameters =
        FDRS.sorted().flatMap { fdr ->
            GAPS.sorted().map { gap -> fdr to gap }
        }

    override val transform: (Pair<Double, Int>) -> String = { (fdr, gap) -> "${fdr}_${gap}" }

    override fun callPeaks(configuration: DataConfig, p: Path, parameter: Pair<Double, Int>) {
        throw IllegalStateException("Batch folder peak calling not valid for $id")
    }

    override fun fileName(cell: Cell, replicate: String, target: String, parameter: Pair<Double, Int>): String {
        return "${cell}_${replicate}_${target}_${SPAN_DEFAULT_BIN}_${parameter.first}_${parameter.second}_peaks.bed"
    }

    override fun defaultParams(uli: Boolean) = SPAN_DEFAULT_FDR to SPAN_DEFAULT_GAP

    override fun tune(
        configuration: DataConfig,
        path: Path,
        target: String,
        useInput: Boolean,
        saveAllPeaks: Boolean
    ) {
        if (!checkTuningRequired(configuration, path, target, useInput)) {
            return
        }
        val folder = folder(path, target, useInput)
        val results = TuningResults()
        val inputPath = if (useInput) configuration.tracksMap.entries
            .filter { it.key.dataType == "chip-seq" && ChipSeqTarget.isInput(it.key.dataType) }
            .flatMap { it.value }.map { it.second.path }.firstOrNull() else null
        val labelledTracks = configuration.extractLabelledTracks(target)

        labelledTracks.forEach { (cellId, replicate, trackPath, labelsPath) ->
            val peakCallingExperiment = SpanPeakCallingExperiment.getExperiment(
                configuration.genomeQuery,
                listOf(SpanDataPaths(trackPath, inputPath!!)),
                SPAN_DEFAULT_BIN, AutoFragment
            )

            if (saveAllPeaks) {
                val progress = Progress {
                    title = "Processing $id $target $cellId $replicate optimal peaks"
                }.bounded(parameters.size.toLong())
                parameters.forEach { parameter ->
                    val peaksPath = folder / transform(parameter) / fileName(cellId, replicate, target, parameter)
                    peaksPath.checkOrRecalculate(ignoreEmptyFile = true) { (path) ->
                        Peak.savePeaks(
                            peakCallingExperiment.results.getFdrGapPeaks(
                                configuration.genomeQuery,
                                parameter.first, parameter.second
                            ), path
                        )
                    }
                    progress.report()
                }
                progress.done()
            }

            LOG.info("Tuning $id $target $cellId $replicate peaks")
            val labels = LocationLabel.loadLabels(
                labelsPath, configuration.genomeQuery.genome
            )
            val (labelErrorsGrid, index) =
                tune(
                    peakCallingExperiment.results, labels,
                    "$id $target $cellId $replicate", parameters,
                    PeaksType.PEAKS_TYPE_ISLANDS
                )


            LOG.info("Saving $id $target $cellId $replicate optimal peaks to $folder")
            cleanup(folder, cellId, replicate)
            val optimalParameters = parameters[index]
            val optimalPeaksPath = folder / fileName(cellId, replicate, target, optimalParameters)
            Peak.savePeaks(
                peakCallingExperiment.results.getFdrGapPeaks(
                    configuration.genomeQuery,
                    optimalParameters.first, optimalParameters.second
                ),
                optimalPeaksPath
            )

            labelErrorsGrid.forEachIndexed { i, error ->
                results.addRecord(replicate, transform(parameters[i]), error, i == index)
            }
        }

        results.saveTuningErrors(folder / "${target}_${id}_errors.csv")
        results.saveOptimalResults(folder / "${target}_${id}_parameters.csv")
        saveGrid(path, target, useInput)
    }

    fun tune(
        results: SpanFitResults,
        labels: List<LocationLabel>,
        id: String,
        parameters: List<Pair<Double, Int>>,
        peaksType: PeaksType,
        cancellableState: CancellableState = CancellableState.current()
    ): Pair<List<LabelErrors>, Int> {
        MultitaskProgress.addTask(id, parameters.size.toLong())
        val tuneResults =
            tune(results, results.fitInfo.genomeQuery(), labels, id, parameters, peaksType, cancellableState)
        MultitaskProgress.finishTask(id)
        return tuneResults
    }


    fun tune(
        results: SpanFitResults,
        genomeQuery: GenomeQuery,
        labels: List<LocationLabel>,
        id: String,
        parameters: List<Pair<Double, Int>>,
        peaksType: PeaksType,
        cancellableState: CancellableState = CancellableState.current()
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
        tuningExecutor.awaitAll(
            parameters.mapIndexed { index, (fdr, gap) ->
                Callable {
                    cancellableState.checkCanceled()
                    val peaksOnLabeledGenomeQuery = results.getPeaks(labeledGenomeQuery, fdr, gap, peaksType)
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
