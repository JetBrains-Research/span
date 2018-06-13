package org.jetbrains.bio.experiments.tuning

import com.google.common.collect.ComparisonChain
import com.google.common.primitives.Shorts
import kotlinx.support.jdk7.use
import org.jetbrains.bio.datasets.CellId
import org.jetbrains.bio.datasets.ChipSeqTarget
import org.jetbrains.bio.datasets.DataType
import org.jetbrains.bio.datasets.toDataType
import org.jetbrains.bio.experiments.DataConfig
import org.jetbrains.bio.experiments.fit.CoverageFitExperiment
import org.jetbrains.bio.experiments.fit.CoverageFitResults
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.span.Span
import org.jetbrains.bio.span.getPeaks
import org.jetbrains.bio.span.savePeaks
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.util.*
import java.nio.file.Path

object SPAN : Tool2Tune<Pair<Double, Int>>() {

    override val id = "span"
    override val suffix = "_peaks.bed"
    const val DEFAULT_BIN = 200
    const val DEFAULT_FDR = 1E-6

    val FDRS = doubleArrayOf(1E-1, 1E-2, 1E-4, DEFAULT_FDR, 1E-8, 1E-10, 1E-12)
    const val DEFAULT_GAP = 5
    val GAPS = intArrayOf(2, DEFAULT_GAP, 10, 20, 40, 60, 80, 100, 120)
    override val parameters =
            FDRS.flatMap { fdr ->
                GAPS.map { gap -> fdr to gap }
            }
    override val transform: (Pair<Double, Int>) -> String = { (fdr, gap) -> "${fdr}_${gap}" }

    override fun callPeaks(configuration: DataConfig, p: Path, parameter: Pair<Double, Int>) {
        throw IllegalStateException("Batch folder peak calling not valid for $id")
    }

    override fun fileName(cellId: CellId, replicate: String, target: String, parameter: Pair<Double, Int>): String {
        return "${cellId}_${replicate}_${target}_${DEFAULT_BIN}_${parameter.first}_${parameter.second}_peaks.bed"
    }

    override fun defaultParams(uli: Boolean) = DEFAULT_FDR to DEFAULT_GAP

    override fun tune(configuration: DataConfig,
                      path: Path,
                      target: String,
                      useInput: Boolean,
                      saveAllPeaks: Boolean) {
        if (!checkTuningRequired(configuration, path, target, useInput)) {
            return
        }
        val folder = folder(path, target, useInput)
        val results = TuningResults()
        val inputPath = if (useInput) configuration.tracksMap.entries
                .filter { it.key.dataType.toDataType() == DataType.CHIP_SEQ && ChipSeqTarget.isInput(it.key.dataType) }
                .flatMap { it.value }.map { it.second.path }.firstOrNull() else null
        val labelledTracks = configuration.extractLabelledTracks(target)

        labelledTracks.forEach { (cellId, replicate, signalPath, labelsPath) ->
            val peakCallingExperiment = Span.getPeakCallingExperiment(
                    configuration.genomeQuery,
                    queries = listOf(ReadsQuery(configuration.genomeQuery, signalPath, inputPath)),
                    bin = SPAN.DEFAULT_BIN)

            if (saveAllPeaks) {
                val progress = Progress {
                    title = "Processing ${SPAN.id} $target $cellId $replicate optimal peaks"
                }.bounded(parameters.size.toLong())
                parameters.forEach { parameter ->
                    val peaksPath = folder / transform(parameter) / fileName(cellId, replicate, target, parameter)
                    peaksPath.checkOrRecalculate(ignoreEmptyFile = true) { (path) ->
                        savePeaks(peakCallingExperiment.getPeaks(parameter.first, parameter.second), path)
                    }
                    progress.report()
                }
                progress.done()
            }

            PeakCallerTuning.LOG.info("Tuning ${SPAN.id} $target $cellId $replicate peaks")
            val labels = PeakAnnotation.loadLabels(labelsPath, configuration.genomeQuery.build)
            val (labelErrorsGrid, index) =
                    tune(peakCallingExperiment, labels, "${SPAN.id} $target $cellId $replicate", parameters)


            PeakCallerTuning.LOG.info("Saving ${SPAN.id} $target $cellId $replicate optimal peaks to $folder")
            cleanup(folder, cellId, replicate)
            val optimalParameters = parameters[index]
            val optimalPeaksPath = folder / fileName(cellId, replicate, target, optimalParameters)
            savePeaks(peakCallingExperiment.getPeaks(optimalParameters.first, optimalParameters.second), optimalPeaksPath)
            labelErrorsGrid.forEachIndexed { i, error ->
                results.addRecord(replicate, transform(parameters[i]), error, i == index)
            }
        }

        results.saveTuningErrors(folder / "${target}_${SPAN}_errors.csv")
        results.saveOptimalResults(folder / "${target}_${SPAN}_parameters.csv")
        saveGrid(path, target, useInput)
    }

    fun <Model : ClassificationModel, State : Any>
            tune(experiment: CoverageFitExperiment<Model, State>,
                 labels: List<PeakAnnotation>,
                 id: String,
                 parameters: List<Pair<Double, Int>>,
                 cancellableState: CancellableState = CancellableState.current()): Pair<List<LabelErrors>, Int> {
        val results = experiment.results
        MultitaskProgress.addTask(id, parameters.size.toLong())
        val tuneResults = tune(results, experiment.genomeQuery, labels, id, parameters, cancellableState)
        MultitaskProgress.finishTask(id)
        return tuneResults
    }


    fun tune(results: CoverageFitResults,
             genomeQuery: GenomeQuery,
             labels: List<PeakAnnotation>,
             id: String,
             parameters: List<Pair<Double, Int>>,
             cancellableState: CancellableState = CancellableState.current())
            : Pair<List<LabelErrors>, Int> {
        val labeledGenomeQuery = genomeQuery.only(labels.map { it.location.chromosome.name }.distinct())
        val labelErrorsGrid = parameters.map { (fdr, gap) ->
            cancellableState.checkCanceled()
            val peaksOnLabeledGenomeQuery = results.getPeaks(labeledGenomeQuery, fdr, gap)
            val errors = computeErrors(labels,
                    LocationsMergingList.create(labeledGenomeQuery, peaksOnLabeledGenomeQuery.map { it.location }))
            MultitaskProgress.reportTask(id)
            errors
        }
        check(labelErrorsGrid.size == parameters.size) {
            "Discrepancy in errors size and parameters size: $labelErrorsGrid vs $parameters"
        }
        val minTotalError = labelErrorsGrid.map { it.error() }.min()!!
        return labelErrorsGrid to parameters.indices
                .filter { labelErrorsGrid[it].error() == minTotalError }
                .sortedWith(Comparator { i1, i2 ->
                    // In case of tie errors prefer smallest FDR and biggest GAP
                    ComparisonChain.start()
                            .compare(parameters[i1].first, parameters[i2].first)
                            .compare(parameters[i2].second, parameters[i1].second).result()
                }).first()
    }

}

object SPAN_REPLICATED : ReplicatedTool2Tune<Pair<Double, Int>>() {

    override val id = "span_replicated"
    override val suffix = SPAN.suffix
    const val DEFAULT_BIN = SPAN.DEFAULT_BIN
    val DEFAULT_FDR = 1E-60
    val FDRS = doubleArrayOf(DEFAULT_FDR, 1E-80, 1E-100, 1E-120, 1E-140, 1E-160)
    val GAPS = SPAN.GAPS
    override val transform = SPAN.transform
    override val parameters =
            FDRS.flatMap { fdr ->
                GAPS.map { gap -> fdr to gap }
            }

    override fun defaultParams(uli: Boolean) = DEFAULT_FDR to SPAN.DEFAULT_GAP

    fun fileName(target: String, parameter: Pair<Double, Int>) =
            "${target}_${DEFAULT_BIN}_${parameter.first}_${parameter.second}_peaks.bed"

    override fun findReplicatedPeaks(path: Path, target: String, useInput: Boolean): Path {
        val folder = folder(path, target, useInput)
        check(folder.exists) {
            "Folder $target $id $folder doesn't exist"
        }
        val peaks = folder.glob("*${target}_*$suffix").firstOrNull()
        check(peaks != null) {
            "No replicated peak file was found for $id and $target at $folder"
        }
        return peaks!!
    }


    /**
     * @return true in case replicate peaks missing of [parameters] has changed
     */
    private fun checkTuningReplicatedRequired(path: Path, target: String, useInput: Boolean): Boolean {
        val grid = parameters.joinToString(separator = ",", transform = transform)
        val existingPeaks = try {
            findReplicatedPeaks(path, target, useInput)
        } catch (ignore: Exception) {
            null
        }
        if (existingPeaks == null || !existingPeaks.exists) {
            PeakCallerTuning.LOG.info("$target $id Missing peak files detected, launching tuning.")
            return true
        }
        if (loadGrid(path, target, useInput) != grid) {
            PeakCallerTuning.LOG.info("$target $id Grid change detected, launching tuning.")
            return true
        }
        PeakCallerTuning.LOG.info("$target $id No action needed, exiting.")
        return false
    }

    override fun tune(configuration: DataConfig,
                      path: Path,
                      target: String,
                      useInput: Boolean,
                      saveAllPeaks: Boolean) {

        if (!checkTuningReplicatedRequired(path, target, useInput)) {
            return
        }
        val folder = folder(path, target, useInput)

        // Get all the replicated tracks
        val labelledTracks = configuration.extractLabelledTracks(target)
        val inputPath = configuration.tracksMap.entries
                .filter { it.key.dataType.toDataType() == DataType.CHIP_SEQ && ChipSeqTarget.isInput(it.key.dataType) }
                .flatMap { it.value }.map { it.second.path }.firstOrNull()

        val replicatedPeakCallingExperiment = Span.getPeakCallingExperiment(
                configuration.genomeQuery,
                queries = labelledTracks.map(LabelledTrack::signalPath).map {
                    ReadsQuery(configuration.genomeQuery, it, inputPath)
                },
                bin = DEFAULT_BIN)

        if (saveAllPeaks) {
            PeakCallerTuning.LOG.info("Saving all $id peaks $target")
            val progress = Progress {
                title = "Processing $id $target"
            }.bounded(parameters.size.toLong())
            SPAN.parameters.forEach { parameter ->
                val peaksPath = folder / transform(parameter) / fileName(target, parameter)
                peaksPath.checkOrRecalculate(ignoreEmptyFile = true) { (path) ->
                    savePeaks(replicatedPeakCallingExperiment.getPeaks(parameter.first, parameter.second), path)
                }
                progress.report()
            }
            progress.done()
        }

        PeakCallerTuning.LOG.info("Tuning $id peaks $target")
        val labelsPath = labelledTracks.first().labelPath
        val labels = PeakAnnotation.loadLabels(labelsPath, configuration.genome)
        val (labelErrorsGrid, index) = SPAN.tune(replicatedPeakCallingExperiment, labels, target, parameters)

        PeakCallerTuning.LOG.info("Saving $target optimal $id peaks to $folder")
        val replicate = "all"
        val existingPaths = findReplicatedPeaks(folder, target, useInput)
        existingPaths.forEach {
            PeakCallerTuning.LOG.info("Removing obsolete file $it")
            it.deleteIfExists()
        }
        savePeaks(replicatedPeakCallingExperiment.getPeaks(parameters[index].first, parameters[index].second),
                folder / fileName(target, parameters[index]))

        val results = TuningResults()
        labelErrorsGrid.forEachIndexed { i, error ->
            results.addRecord(replicate, transform(parameters[i]), error, index == i)
        }
        results.saveTuningErrors(folder / "${target}_${SPAN_REPLICATED}_errors.csv")
        results.saveOptimalResults(folder / "${target}_${SPAN_REPLICATED}_parameters.csv")

        val labelErrorsPath = folder / "${target}_${SPAN_REPLICATED}_error_labels.bed"
        BedFormat().print(labelErrorsPath).use { printer ->
            for (entry in labelErrorsGrid[index]) {
                val score = Shorts.checkedCast((entry.value.rate() * 1000).toLong())
                printer.print(entry.key.asBedEntry().copy(score = score))
            }
        }
        saveGrid(path, target, useInput)
    }
}