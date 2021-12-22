package org.jetbrains.bio.experiments.tuning

import kotlinx.support.jdk7.use
import org.jetbrains.bio.experiments.tuning.tools.Span
import org.jetbrains.bio.genome.data.Cell
import org.jetbrains.bio.genome.data.DataConfig
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.nio.file.StandardOpenOption


abstract class Tool2Tune<T> {

    abstract val id: String

    abstract val suffix: String

    /**
     * In case of equal tuning error, first parameter will be chosen,
     * so that each tool has its own order of parameters.
     */
    abstract val parameters: List<T>

    abstract val transform: (T) -> String

    abstract fun callPeaks(configuration: DataConfig, p: Path, parameter: T)

    abstract fun fileName(cell: Cell, replicate: String, target: String, parameter: T): String

    abstract fun defaultParams(uli: Boolean): T

    override fun toString(): String {
        return id
    }

    /**
     * Returns folder for given [target] and [tool], takes into account [useInput] in case of [Span]
     */
    fun folder(path: Path, target: String, useInput: Boolean) =
        (path / target /
                (if (this == Span)
                    "$id${if (!useInput) "_noinput" else ""}"
                else
                    id)).apply {
            createDirectories()
        }

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

    protected fun saveGrid(path: Path, target: String, useInput: Boolean) {
        gridFile(path, target, useInput)
            .bufferedWriter(StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING)
            .use { it.write(parameters.joinToString(separator = ",", transform = transform)) }
    }

    /**
     * @return map of Pair<cell, replicate> -> path
     */
    fun tunedPeaks(
        configuration: DataConfig,
        path: Path,
        target: String,
        useInput: Boolean
    ): Map<Pair<Cell, String>, Path> {
        val folder = folder(path, target, useInput)
        if (folder.notExists) {
            LOG.warn("Folder $target $this $folder doesn't exist")

            return emptyMap()
        }
        val labelledTracks = configuration.extractLabelledTracks(target)
        val paths = labelledTracks.map {
            (it.cell to it.name) to
                    findPeakFiles(folder, it.cell.toString(), it.name).firstOrNull()
        }
        val existingPaths = paths.filter { it.second != null }.map { (a, b) -> a to b!! }
        if (existingPaths.size != paths.size) {
            LOG.warn("Not all peak files were found for $id and $target at $folder\n" +
                    "Missing paths for : ${paths.filter { it.second == null }}"
            )
        }
        return existingPaths.associate { it }
    }


    /**
     * @return true in case when some peaks are missing of [parameters] has changed
     */
    protected fun checkTuningRequired(
        configuration: DataConfig,
        path: Path,
        target: String,
        useInput: Boolean
    ): Boolean {
        val prefix = "$target $id"
        val existingPeaks = tunedPeaks(configuration, path, target, useInput)
        val labelledTracks = configuration.extractLabelledTracks(target)
        if (labelledTracks.size != existingPeaks.size) {
            LOG.info("$prefix Missing peak files detected, launching tuning.")
            return true
        }
        if (loadGrid(path, target, useInput) != parameters.joinToString(separator = ",", transform = transform)) {
            LOG.info("$prefix Grid change detected, launching tuning.")
            return true
        }
        LOG.info("$prefix No action needed, exiting.")
        return false
    }

    abstract fun tune(
        configuration: DataConfig,
        path: Path,
        target: String,
        useInput: Boolean,
        saveAllPeaks: Boolean
    )


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
            LOG.warn("NO peaks folder for $folder $id $cellId $donor")
            emptyList()
        }
    }

    protected fun cleanup(folder: Path, cellId: Cell, replicate: String) {
        findPeakFiles(folder, cellId.name, replicate).forEach {
            LOG.info("Removing obsolete file $it")
            it.deleteIfExists()
            "${it}_rip.csv".toPath().deleteIfExists()
        }
    }

    companion object {
        private val LOG: Logger = LoggerFactory.getLogger(Tool2Tune::class.java)
    }
}

abstract class ReplicatedTool2Tune<T> : Tool2Tune<T>() {

    override fun callPeaks(configuration: DataConfig, p: Path, parameter: T) {
        throw IllegalStateException("Batch #callPeaks is not available in ${id}")
    }

    override fun fileName(cell: Cell, replicate: String, target: String, parameter: T): String {
        throw IllegalStateException("#fileName is not available in ${id}")
    }

    abstract fun findReplicatedPeaks(path: Path, target: String, useInput: Boolean): Path

}