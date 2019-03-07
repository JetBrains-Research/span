package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.dataset.ChipSeqTarget
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.util.checkOrRecalculate
import org.jetbrains.bio.util.delete
import org.jetbrains.bio.util.div
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.util.*

/**
 * Creates symbolic links in the [basePath] given [modification]
 */
fun DataConfig.fillData(
        basePath: Path,
        modification: String?,
        addInput: Boolean = true,
        addFailedTracks: Boolean = false
): List<Path> {
    val files = ArrayList<Path>()
    for ((key, section) in tracksMap) {
        for ((replicate, contents) in section.filter { addFailedTracks || !it.second.failedTrack }) {
            val isInput = ChipSeqTarget.isInput(key.dataType)
            val addTrack = if (isInput) addInput else (modification == null || modification == key.dataType)
            if (addTrack) {
                val localName = "${key.cellId}_${
                if (isInput)
                    "input"
                else
                    "${replicate}_${key.dataType}"
                }.bam"
                val localPath = basePath / localName
                if (Files.exists(localPath, LinkOption.NOFOLLOW_LINKS)) {
                    localPath.delete()
                }
                DataConfig.LOG.info("${contents.path} -> $localPath")
                Files.createSymbolicLink(localPath, contents.path)
                files.add(localPath)
            }
        }
    }
    return files
}

/**
 * Creates symbolic links using [fillData] and launches [body] code.
 */
fun DataConfig.runBatch(
        path: Path, mark: String?,
        addInput: Boolean = true,
        addFailedTracks: Boolean = false,
        body: (Path) -> Unit
) {
    path.checkOrRecalculate(isDirectory = true) {
        val basePath = it.path
        val files = fillData(basePath, mark, addInput = addInput, addFailedTracks = addFailedTracks)
        body(basePath)
        files.forEach { Files.delete(it) }
    }
}