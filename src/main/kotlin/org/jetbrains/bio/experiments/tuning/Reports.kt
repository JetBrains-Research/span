package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.coverage.removeDuplicates
import org.jetbrains.bio.dataframe.DataFrameBuilder
import org.jetbrains.bio.dataframe.DataFrameSpec
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.tools.Washu
import org.jetbrains.bio.util.*
import java.nio.file.Path

/**
 * Methods to create summary peak caller tuning data frame
 */

internal fun report(): DataFrameBuilder {
    return DataFrameBuilder(DataFrameSpec().strings("donor").strings("modification").strings("tool")
            .ints("peaks").ints("length").floats("frip")
            .strings("procedure").strings("params").strings("file").strings("status"))
}


private operator fun <T> List<T>.component6(): T = get(5)


private fun findPeaks(peaksPath: Path): List<Path> {
    return (peaksPath.glob("*_peaks.bed") +
            peaksPath.glob("*island.bed") +
            peaksPath.glob("*.broadPeak") +
            peaksPath.glob("*.narrowPeak"))
}

private fun computeFRIP(configuration: DataConfig, target: String, folder: Path, washu: Washu) {
    val folderPeaks = findPeaks(folder)

    configuration.tracksMap.entries
            .filter { it.key.dataType == target }
            .flatMap { entry -> entry.value }.forEach { (_, replicateData) ->
                val uniqueBamPath = removeDuplicates(replicateData.path)
                val donor = "GSM\\d+|[YO]D\\d+".toRegex().find(replicateData.path.toString())!!.value
                val peaks = folderPeaks.filter { "${donor}_" in it.stem }
                check(peaks.isNotEmpty()) {
                    "No peak files found for $donor in $folder"
                }
                check(peaks.size == 1) {
                    Washu.LOG.warn("More than 1 peak file found for $donor in $folder\n$peaks")
                }
                washu.runRip(uniqueBamPath, peaks.first())
            }
}

internal fun PeakCallerTuning.computeFripAndReport(report: DataFrameBuilder,
                                                   target: String,
                                                   tool: Tool2Tune<*>,
                                                   folder: Path,
                                                   procedure: String,
                                                   washu: Washu) {
    computeFRIP(configuration, target, folder, washu)

    val toolName = tool.id
    val failed = configuration.tracksMap.filter { it.key.dataType == target }
            .flatMap { it.value.filter { it.second.failedTrack }.map { it.first } }.toSet()
    findPeaks(folder).forEach { path ->
        PeakCallerTuning.LOG.info("$target $toolName $procedure $path")
        val donor = donor(path)
        val params = params(path, target)
        val peaks = path.useLines { it.count() }
        var length = 0
        var frip = 0f
        val ripCsv = (path.toString() + "_rip.csv").toPath()
        if (ripCsv.exists) {
            val ripLine = ripCsv.bufferedReader().useLines {
                it.last()
            }
            try {
                val (_, _, reads, ps, l, rip) = ripLine.split(",")
                check(peaks == ps.toInt()) {
                    "Different peaks number in $path and $ripCsv"
                }
                if (peaks > 0) {
                    frip = rip.toFloat() / reads.toLong()
                    length = l.toInt()
                }
            } catch (e: Throwable) {
                PeakCallerTuning.LOG.error("Failed to analyse $ripCsv\n$ripLine", e)
            }
        } else {
            path.useLines {
                it.forEach {
                    val split = it.split("\t")
                    length += split[2].toInt() - split[1].toInt()
                }
            }
        }
        report.add(donor, target, toolName,
                peaks, length, frip,
                procedure, params, path.toAbsolutePath().toString(),
                if (donor in failed) "failed" else "ok")
    }
}

/**
 * See [DataConfig.fillData] for details
 */
internal fun donor(path: Path) = path.stem.split('_')[1]

/**
 * See [DataConfig.fillData] for details on naming
 */
internal fun params(path: Path, target: String) =
        path.stem.substringAfter(target).replace("peaks|island".toRegex(), "").trim('_', '-')


