package org.jetbrains.bio.experiments.tuning

import kotlinx.support.jdk7.use
import org.jetbrains.bio.experiments.tuning.PeakCallerTuning.Companion.LOG
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.util.*
import java.nio.file.Path
import kotlin.math.roundToInt

/**
 * Created by Aleksei Dievskii on 15.02.2018.
 */

internal fun PeakCallerTuning.generateLabelErrors(target: String, tool: Tool2Tune<*>): Path? {
    val labelErrorsPath = tool.folder(experimentPath, target, useInput) / "${target}_${tool}_error_labels.bed"
    labelErrorsPath.checkOrRecalculate(ignoreEmptyFile = true) { output ->
        val tunedPeaks = tool.tunedPeaks(configuration, experimentPath, target, useInput)
        val labelledTracks = configuration.extractLabelledTracks(target)
        val labelErrors = LabelErrors()
        labelledTracks.forEach { track ->
            val labels = PeakAnnotation.loadLabels(track.labelPath, configuration.genomeQuery.genome)
            val peaksPath = tunedPeaks[track.cellId to track.name]
            if (peaksPath == null) {
                LOG.warn("No tuned peaks for $target, $tool, ${track.cellId} and ${track.name}, " +
                        "can't generate label errors")
                return@checkOrRecalculate
            }
            labelErrors.combine(computeErrors(labels, LocationsMergingList.load(configuration.genomeQuery, peaksPath)))
        }
        BedFormat().print(output.path).use { printer ->
            labelErrors.forEach { entry ->
                val score = (entry.value.rate() * 1000).roundToInt()
                printer.print(entry.key.asBedEntry().copy(score = score))
            }
        }
    }
    return if (labelErrorsPath.exists && labelErrorsPath.size.bytes > 0) {
        labelErrorsPath
    } else {
        labelErrorsPath.deleteIfExists()
        null
    }
}

internal fun PeakCallerTuning.generateReplicatedLabelErrors(target: String, tool: ReplicatedTool2Tune<*>): Path? {
    val labelErrorsPath = tool.folder(experimentPath, target, useInput) / "${target}_${tool}_error_labels.bed"
    labelErrorsPath.checkOrRecalculate(ignoreEmptyFile = true) { output ->
        val peaksPath = try {
            tool.findReplicatedPeaks(experimentPath, target, useInput)
        } catch (e: IllegalStateException) {
            null
        }
        val labels = PeakAnnotation.loadLabels(configuration.extractLabelledTracks(target).first().labelPath,
                configuration.genomeQuery.genome)
        if (peaksPath == null) {
            LOG.warn("No tuned replicated peak file for $target and $tool, can't generate label errors")
            return@checkOrRecalculate
        }
        val labelErrors = computeErrors(labels, LocationsMergingList.load(configuration.genomeQuery, peaksPath))
        BedFormat().print(output.path).use { printer ->
            labelErrors.forEach { entry ->
                val score = (entry.value.rate() * 1000).roundToInt()
                printer.print(entry.key.asBedEntry().copy(score = score))
            }
        }
    }
    return if (labelErrorsPath.exists && labelErrorsPath.size.bytes > 0) {
        labelErrorsPath
    } else {
        labelErrorsPath.deleteIfExists()
        null
    }
}