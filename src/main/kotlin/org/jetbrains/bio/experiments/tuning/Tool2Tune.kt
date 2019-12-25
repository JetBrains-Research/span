package org.jetbrains.bio.experiments.tuning

import kotlinx.support.jdk7.use
import org.apache.log4j.Logger
import org.jetbrains.bio.dataset.CellId
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.experiments.tuning.tools.Span
import org.jetbrains.bio.util.*
import java.nio.file.Path
import java.nio.file.StandardOpenOption


abstract class Tool2Tune<T> {

    private val LOG = Logger.getLogger(Tool2Tune::class.java)

    abstract val id: String

    abstract val suffix: String

    /**
     * In case of equal tuning error, first parameter will be chosen,
     * so that each tool has it's own order of parameters.
     */
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
    fun tunedPeaks(configuration: DataConfig,
                   path: Path,
                   target: String,
                   useInput: Boolean): Map<Pair<CellId, String>, Path> {
        val folder = folder(path, target, useInput)
        if (folder.notExists) {
            LOG.warn("Folder $target $this $folder doesn't exist")

            return emptyMap()
        }
        val labelledTracks = configuration.extractLabelledTracks(target)
        val paths = labelledTracks.map {
            (it.cellId to it.name) to
                    findPeakFiles(folder, it.cellId.toString(), it.name).firstOrNull()
        }
        val existingPaths = paths.filter { it.second != null }.map { (a, b) -> a to b!! }
        if (existingPaths.size != paths.size) {
            LOG.warn("Not all peak files were found for $id and $target at $folder\n" +
                    "Missing paths for : ${paths.filter { it.second == null }}")
        }
        return existingPaths.associate { it }
    }


    /**
     * @return true in case when some peaks are missing of [parameters] has changed
     */
    protected fun checkTuningRequired(configuration: DataConfig,
                                      path: Path,
                                      target: String,
                                      useInput: Boolean): Boolean {
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
            saveAllPeaks: Boolean)


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

    protected fun cleanup(folder: Path, cellId: CellId, replicate: String) {
        findPeakFiles(folder, cellId.name, replicate).forEach {
            LOG.info("Removing obsolete file $it")
            it.deleteIfExists()
            "${it}_rip.csv".toPath().deleteIfExists()
        }
    }

    companion object {
        /**
         * Combines all the combinations of parameters starting from defaults in each params
         * @param params: pairs of list of params, default value index.
         */
        fun combineParams(vararg params: Pair<List<Any>, Int>): List<Array<Any>> {
            if (params.isEmpty()) {
                return emptyList()
            }
            val result = arrayListOf<Array<Any>>()
            val tuple = Array<Any>(params.size) { 0 }
            processIndex(params, 0, tuple, result)
            return result
        }

        private fun processIndex(params: Array<out Pair<List<Any>, Int>>,
                                 index: Int,
                                 tuple: Array<Any>,
                                 result: ArrayList<Array<Any>>) {
            val (paramsI, defaultIndexI) = params[index]
            var indexOffset = 0
            while (defaultIndexI + indexOffset in paramsI.indices || defaultIndexI - indexOffset in paramsI.indices) {
                val pIndex = defaultIndexI + indexOffset
                if (pIndex in paramsI.indices) {
                    tuple[index] = paramsI[pIndex]
                    if (index == params.size - 1) {
                        result.add(tuple.clone())
                    } else {
                        processIndex(params, index + 1, tuple, result)
                    }
                }
                indexOffset = (if (indexOffset <= 0) -indexOffset + 1 else -indexOffset)
            }
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