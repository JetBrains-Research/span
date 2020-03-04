package org.jetbrains.bio.experiments.tuning.tools

import kotlinx.support.jdk7.use
import org.apache.log4j.Logger
import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.dataset.ChipSeqTarget
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.dataset.DataType
import org.jetbrains.bio.dataset.toDataType
import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_BIN
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_GAP
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment.Companion.SPAN_REPLICATED_DEFAULT_FDR
import org.jetbrains.bio.experiments.tuning.*
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.span.getPeaks
import org.jetbrains.bio.span.savePeaks
import org.jetbrains.bio.util.*
import java.nio.file.Path
import kotlin.math.roundToInt

object SpanReplicated : ReplicatedTool2Tune<Pair<Double, Int>>() {

    private val LOG = Logger.getLogger(SpanReplicated::class.java)

    override val id = "span_replicated"
    override val suffix = Span.suffix

    private val FDRS = listOf(0.05, 0.01, 1E-4, 1E-6, 1E-8, SPAN_REPLICATED_DEFAULT_FDR,
            1E-20, 1E-40, 1E-60, 1E-80, 1E-100)

    private val GAPS = Span.GAPS

    override val transform = Span.transform

    override val parameters =
            FDRS.sorted().flatMap { fdr ->
                GAPS.sorted().map { gap -> fdr to gap }
            }

    override fun defaultParams(uli: Boolean) = SPAN_REPLICATED_DEFAULT_FDR to SPAN_DEFAULT_GAP

    fun fileName(target: String, parameter: Pair<Double, Int>) =
            "${target}_${SPAN_DEFAULT_BIN}_${parameter.first}_${parameter.second}_peaks.bed"

    override fun findReplicatedPeaks(path: Path, target: String, useInput: Boolean): Path {
        val peaks = findReplicatedFiles(path, target, useInput)
        check(peaks.size == 1) {
            "No or more than 1 replicated peaks file was found for $id and $target at $path"
        }
        return peaks.first()
    }

    private fun findReplicatedFiles(path: Path, target: String, useInput: Boolean): List<Path> {
        val folder = folder(path, target, useInput)
        check(folder.exists) {
            "Folder $target $id $folder doesn't exist"
        }
        return folder.glob("*${target}_*$suffix").toList()
    }


    /**
     * @return true in case replicate peaks missing of [parameters] has changed
     */
    private fun checkTuningReplicatedRequired(path: Path, target: String, useInput: Boolean): Boolean {
        val grid = parameters.joinToString(separator = ",", transform = transform)
        val existingPeaks = findReplicatedFiles(path, target, useInput)
        if (existingPeaks.size != 1) {
            LOG.info("$target $id Missing or extra peak files detected, launching tuning.")
            existingPeaks.forEach {
                LOG.info("Removing obsolete file $it")
                it.deleteIfExists()
            }
            return true
        }
        if (loadGrid(path, target, useInput) != grid) {
            LOG.info("$target $id Grid change detected, launching tuning.")
            return true
        }
        LOG.info("$target $id No action needed, exiting.")
        return false
    }

    override fun tune(
            configuration: DataConfig,
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

        val replicatedPeakCallingExperiment = SpanPeakCallingExperiment.getExperiment(
                configuration.genomeQuery,
                labelledTracks.map(LabelledTrack::trackPath).map { SpanDataPaths(it, inputPath) },
                SPAN_DEFAULT_BIN, AutoFragment
        )

        if (saveAllPeaks) {
            LOG.info("Saving all $id peaks $target")
            val progress = Progress {
                title = "Processing $id $target"
            }.bounded(parameters.size.toLong())
            parameters.forEach { parameter ->
                val peaksPath = folder / transform(parameter) / fileName(target, parameter)
                peaksPath.checkOrRecalculate(ignoreEmptyFile = true) { (p) ->
                    savePeaks(replicatedPeakCallingExperiment.results.getPeaks(configuration.genomeQuery,
                            parameter.first, parameter.second),
                            p)
                }
                progress.report()
            }
            progress.done()
        }

        LOG.info("Tuning $id peaks $target")
        val labelsPath = labelledTracks.first().labelPath
        val labels = LocationLabel.loadLabels(labelsPath, configuration.genomeQuery.genome)
        val (labelErrorsGrid, index) = Span.tune(replicatedPeakCallingExperiment.results, labels, target, parameters)

        LOG.info("Saving $target optimal $id peaks to $folder")
        savePeaks(replicatedPeakCallingExperiment.results.getPeaks(configuration.genomeQuery,
                parameters[index].first, parameters[index].second),
                folder / fileName(target, parameters[index]))

        val results = TuningResults()
        labelErrorsGrid.forEachIndexed { i, error ->
            results.addRecord("", transform(parameters[i]), error, index == i)
        }
        results.saveTuningErrors(folder / "${target}_${id}_errors.csv")
        results.saveOptimalResults(folder / "${target}_${id}_parameters.csv")

        val labelErrorsPath = folder / "${target}_${id}_error_labels.bed"
        BedFormat().print(labelErrorsPath).use { printer ->
            for (entry in labelErrorsGrid[index]) {
                val score = (entry.value.rate() * 1000).roundToInt()
                printer.print(entry.key.asBedEntry().copy(score = score))
            }
        }
        saveGrid(path, target, useInput)
    }
}

