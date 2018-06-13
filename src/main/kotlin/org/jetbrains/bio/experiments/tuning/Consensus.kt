package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.datasets.CellId
import org.jetbrains.bio.experiments.DataConfig
import org.jetbrains.bio.util.checkOrRecalculate
import org.jetbrains.bio.util.createDirectories
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.exists
import java.nio.file.Path

fun DataConfig.medianConsensus(target: String, cell: CellId? = null, tool: Tool2Tune<*>? = null): List<Path> =
        (if (tool != null) listOf(tool) else INDIVIDUAL_TOOLS).mapNotNull { t ->
            PeakCallerTuning(this, tools = arrayListOf(t)).medianConsensus(target, cell, t)
        }

fun DataConfig.weakConsensus(target: String, cell: CellId? = null, tool: Tool2Tune<*>? = null): List<Path> =
        (if (tool != null) listOf(tool) else INDIVIDUAL_TOOLS).mapNotNull { t ->
            PeakCallerTuning(this, tools = arrayListOf(t)).weakConsensus(target, cell, t)
        }

internal fun PeakCallerTuning.generateConsensus(target: String, tool: Tool2Tune<*>) {
    medianConsensus(target, tool = tool)
    weakConsensus(target, tool = tool)
    configuration.extractLabelledTracks(target).map { it.cellId }.toSet().forEach {
        medianConsensus(target, it, tool)
        weakConsensus(target, it, tool)
    }
}


private fun PeakCallerTuning.consensusFolder(target: String, tool: Tool2Tune<*>) =
        tool.folder(experimentPath, target, useInput) / "consensus"

private fun PeakCallerTuning.medianConsensus(target: String, cell: CellId? = null, tool: Tool2Tune<*>): Path? {
    val consensusFolder = consensusFolder(target, tool)
    consensusFolder.createDirectories()
    val consensusFile = consensusFolder /
            "${target}_${tool}${if (cell != null) "_${cell.name}" else ""}_median_consensus.bed"
    val goodTracks = configuration.extractLabelledTracks(target)
            .filter { !it.failed && (cell == null || it.cellId == cell) }.map { it.cellId to it.name }
    val peaks = tool.tunedPeaks(configuration, experimentPath, target, useInput).filter { it.key in goodTracks }.map { it.value }
    if (!consensusFile.exists && goodTracks.size != peaks.size) {
        PeakCallerTuning.LOG.warn("Couldn't find or generate median consensus for $target $tool, missing peak files.")
        return null
    }
    consensusFile.checkOrRecalculate(ignoreEmptyFile = true) { output ->
        washu.medianNucleotideConsensus(output.path, peaks)
    }
    return consensusFile
}

private fun PeakCallerTuning.weakConsensus(target: String, cell: CellId? = null, tool: Tool2Tune<*>): Path? {
    val consensusFolder = consensusFolder(target, tool)
    consensusFolder.createDirectories()
    val consensusFile = consensusFolder /
            "${target}_${tool}${if (cell != null) "_${cell.name}" else ""}_weak_consensus.bed"
    val goodTracks = configuration.extractLabelledTracks(target)
            .filter { !it.failed && (cell == null || it.cellId == cell) }.map { it.cellId to it.name }
    val peaks = tool.tunedPeaks(configuration, experimentPath, target, useInput).filter { it.key in goodTracks }.map { it.value }
    if (!consensusFile.exists && goodTracks.size != peaks.size) {
        PeakCallerTuning.LOG.warn("Couldn't find or generate weak consensus for $target $tool, missing peak files.")
        return null
    }
    consensusFile.checkOrRecalculate(ignoreEmptyFile = true) { output ->
        washu.weakNucleotideConsensus(output.path, peaks)
    }
    return consensusFile
}
