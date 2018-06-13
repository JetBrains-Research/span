package org.jetbrains.bio.experiments.tuning

import kotlinx.support.jdk7.use
import org.jetbrains.bio.datasets.CellId
import org.jetbrains.bio.experiments.DataConfig
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.span.Span
import org.jetbrains.bio.tools.runBatch
import org.jetbrains.bio.util.*
import java.nio.file.Path
import java.nio.file.StandardCopyOption
import java.nio.file.StandardOpenOption


abstract class Tool2Tune<T> {
    abstract val id: String
    abstract val suffix: String
    abstract val parameters: List<T>
    abstract val transform: (T) -> String

    abstract fun callPeaks(configuration: DataConfig, p: Path, parameter: T)

    abstract fun fileName(cellId: CellId, replicate: String, target: String, parameter: T): String

    abstract fun defaultParams(uli: Boolean): T

    override fun toString(): String {
        return id
    }

    /**
     * Returns folder for given [target] and [tool], takes into account [useInput] in case of [Span]
     */
    internal fun folder(path: Path, target: String, useInput: Boolean) =
            path / target /
                    (if (this == SPAN)
                        "$id${if (!useInput) "_noinput" else ""}"
                    else
                        id)

    fun defaultsFolder(path: Path, target: String, useInput: Boolean, uli: Boolean) =
            folder(path, target, useInput) / transform(defaultParams(uli))

    private fun gridFile(path: Path, target: String, useInput: Boolean) =
            folder(path, target, useInput) / ".grid"

    internal fun loadGrid(path: Path, target: String, useInput: Boolean): String? {
        val gridFile = gridFile(path, target, useInput)
        return if (gridFile.exists && gridFile.isReadable) {
            gridFile.useLines { it.firstOrNull() }
        } else {
            null
        }
    }

    internal fun saveGrid(path: Path, target: String, useInput: Boolean) {
        gridFile(path, target, useInput)
                .bufferedWriter(StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
                .use { it.write(parameters.joinToString(separator = ",", transform = transform)) }
    }

    /**
     * @return map of Pair<cell, replicate> -> path
     */
    fun tunedPeaks(configuration: DataConfig,
                   path: Path,
                   target: String,
                   useInput: Boolean): Map<Pair<CellId, String>, Path> {
        val folder = folder(path, target, useInput)
        if (folder.notExists) {
            PeakCallerTuning.LOG.warn("Folder $target $this $folder doesn't exist")

            return emptyMap()
        }
        val labelledTracks = configuration.extractLabelledTracks(target)
        val paths = labelledTracks.map {
            (it.cellId to it.name) to
                    findPeakFiles(folder, it.cellId.toString(), it.name).firstOrNull()
        }
        val existingPaths = paths.filter { it.second != null }.map { (a, b) -> a to b!! }
        if (existingPaths.size != paths.size) {
            PeakCallerTuning.LOG.warn("Not all peak files were found for $id and $target at $folder\n" +
                    "Missing paths for : ${paths.filter { it.second == null }}")
        }
        return existingPaths.associate { it }
    }


    /**
     * @return true in case when some peaks are missing of [parameters] has changed
     */
    internal fun checkTuningRequired(configuration: DataConfig,
                                     path: Path,
                                     target: String,
                                     useInput: Boolean): Boolean {
        val prefix = "$target $id"
        val existingPeaks = tunedPeaks(configuration, path, target, useInput)
        val labelledTracks = configuration.extractLabelledTracks(target)
        if (labelledTracks.size != existingPeaks.size) {
            PeakCallerTuning.LOG.info("$prefix Missing peak files detected, launching tuning.")
            return true
        }
        if (loadGrid(path, target, useInput) != parameters.joinToString(separator = ",", transform = transform)) {
            PeakCallerTuning.LOG.info("$prefix Grid change detected, launching tuning.")
            return true
        }
        PeakCallerTuning.LOG.info("$prefix No action needed, exiting.")
        return false
    }

    open fun tune(configuration: DataConfig,
                  path: Path,
                  target: String,
                  useInput: Boolean,
                  saveAllPeaks: Boolean) {
        check(saveAllPeaks) {
            "Tuning creates all intermediate peak files by default"
        }
        if (!checkTuningRequired(configuration, path, target, useInput)) {
            return
        }
        val folder = folder(path, target, useInput)

        for (parameter in parameters) {
            PeakCallerTuning.LOG.info("Processing $target $id $parameter")
            configuration.runBatch(folder / transform(parameter), target,
                    addFailedTracks = true,
                    addInput = useInput) { p -> callPeaks(configuration, p, parameter) }
        }

        val results = TuningResults()
        val labelledTracks = configuration.extractLabelledTracks(target)
        val progress = Progress {
            title = "Tuning procedure $target"
        }.bounded(labelledTracks.size.toLong() * parameters.size)
        val labelErrors = LabelErrors()
        for ((cellId, replicate, _, labelsPath) in labelledTracks) {
            val labels = PeakAnnotation.loadLabels(labelsPath, configuration.genome)
            val labelErrorsGrid = parameters.map { parameter ->
                val peaksPath = folder / transform(parameter) /
                        fileName(cellId, replicate, target, parameter)
                val errors = computeErrors(labels, LocationsMergingList.load(configuration.genomeQuery, peaksPath))
                progress.report()
                errors
            }

            val minTotalError = labelErrorsGrid.map { it.error() }.min()!!
            // Take something in the middle to get not that extreme settings
            val optimalIndex = parameters.indices
                    .filter { labelErrorsGrid[it].error() == minTotalError }.let { it[it.size / 2] }
            val optimalParameter = parameters[optimalIndex]
            labelErrors.combine(labelErrorsGrid[optimalIndex])
            cleanup(folder, cellId, replicate)
            // Copy optimal parameters path to folder
            val optimalPeaksPath = folder / transform(optimalParameter) /
                    fileName(cellId, replicate, target, optimalParameter)
            optimalPeaksPath.copy(folder, StandardCopyOption.REPLACE_EXISTING)
            (optimalPeaksPath.toString() + "_rip.csv").toPath().copy(folder, StandardCopyOption.REPLACE_EXISTING)

            labelErrorsGrid.forEachIndexed { i, error ->
                results.addRecord(replicate,
                        transform(parameters[i]),
                        error,
                        minTotalError == error.error())
            }
        }
        progress.done()
        results.saveTuningErrors(folder / "${target}_${id}_errors.csv")
        results.saveOptimalResults(folder / "${target}_${id}_parameters.csv")
        saveGrid(path, target, useInput)
    }


    private fun findPeakFiles(folder: Path, cellId: String, donor: String?): List<Path> {
        return if (folder.exists)
            folder.list().filter {
                it.fileName.toString().let {
                    it.contains(Regex("[^a-zA-Z0-9]$donor[^a-zA-Z0-9]"))
                            && it.contains(cellId)
                            && it.endsWith(suffix)
                }
            }
        else {
            PeakCallerTuning.LOG.warn("NO peaks folder for $folder $id $cellId $donor")
            emptyList()
        }
    }

    internal fun cleanup(folder: Path, cellId: CellId, replicate: String) {
        findPeakFiles(folder, cellId.name, replicate).forEach {
            PeakCallerTuning.LOG.info("Removing obsolete file $it")
            it.deleteIfExists()
            "${it}_rip.csv".toPath().deleteIfExists()
        }
    }

}

abstract class ReplicatedTool2Tune<T> : Tool2Tune<T>() {
    override fun callPeaks(configuration: DataConfig, p: Path, parameter: T) {
        throw IllegalStateException("Batch #callPeaks is not available in ${id}")
    }

    override fun fileName(cellId: CellId, replicate: String, target: String, parameter: T): String {
        throw IllegalStateException("#fileName is not available in ${id}")
    }

    abstract fun findReplicatedPeaks(path: Path, target: String, useInput: Boolean): Path

}

val INDIVIDUAL_TOOLS = listOf(MACS2_BROAD, MACS2, SICER, SPAN)
val REPLICATED_TOOLS = listOf(SPAN_REPLICATED)

