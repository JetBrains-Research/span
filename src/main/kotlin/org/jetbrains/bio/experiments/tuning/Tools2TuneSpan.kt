package org.jetbrains.bio.experiments.tuning

import kotlinx.support.jdk7.use
import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.removeDuplicates
import org.jetbrains.bio.dataset.*
import org.jetbrains.bio.experiments.fit.SpanFitResults
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.span.getPeaks
import org.jetbrains.bio.span.savePeaks
import org.jetbrains.bio.tools.ToolsChipSeqWashu
import org.jetbrains.bio.util.*
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.Executors
import kotlin.math.roundToInt

object Span : Tool2Tune<Pair<Double, Int>>() {

    private val tuningExecutor = Executors.newWorkStealingPool(parallelismLevel())

    override val id = "span"
    override val suffix = "_peaks.bed"
    const val DEFAULT_BIN = 200

    const val DEFAULT_FDR = 1E-6
    private val FDRS = doubleArrayOf(0.1, 0.05, 1E-2, 1E-3, 1E-4, 1E-5, DEFAULT_FDR, 1E-7, 1E-8, 1E-9, 1E-10)

    const val DEFAULT_GAP = 5
    val GAPS = intArrayOf(0, 1, 2, 3, 4, DEFAULT_GAP, 6, 7, 8, 9, 10)

    override val parameters =
            FDRS.sorted().flatMap { fdr ->
                GAPS.sorted().map { gap -> fdr to gap }
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

        labelledTracks.forEach { (cellId, replicate, trackPath, labelsPath) ->
            val peakCallingExperiment = SpanPeakCallingExperiment.getExperiment(configuration.genomeQuery,
                listOf(trackPath to inputPath),
                DEFAULT_BIN, AutoFragment
            )

            if (saveAllPeaks) {
                val progress = Progress {
                    title = "Processing $id $target $cellId $replicate optimal peaks"
                }.bounded(parameters.size.toLong())
                parameters.forEach { parameter ->
                    val peaksPath = folder / transform(parameter) / fileName(cellId, replicate, target, parameter)
                    peaksPath.checkOrRecalculate(ignoreEmptyFile = true) { (path) ->
                        savePeaks(peakCallingExperiment.results.getPeaks(configuration.genomeQuery,
                            parameter.first, parameter.second), path)
                    }
                    progress.report()
                }
                progress.done()
            }

            PeakCallerTuning.LOG.info("Tuning $id $target $cellId $replicate peaks")
            val labels = PeakAnnotation.loadLabels(
                labelsPath, configuration.genomeQuery.genome
            )
            val (labelErrorsGrid, index) =
                    tune(peakCallingExperiment.results, labels, "$id $target $cellId $replicate", parameters)


            PeakCallerTuning.LOG.info("Saving $id $target $cellId $replicate optimal peaks to $folder")
            cleanup(folder, cellId, replicate)
            val optimalParameters = parameters[index]
            val optimalPeaksPath = folder / fileName(cellId, replicate, target, optimalParameters)
            savePeaks(peakCallingExperiment.results.getPeaks(configuration.genomeQuery,
                optimalParameters.first, optimalParameters.second),
                optimalPeaksPath)

            // Compute _rip.sh file
            ToolsChipSeqWashu().runRip(removeDuplicates(trackPath), optimalPeaksPath)
            check("${optimalPeaksPath}_rip.csv".toPath().exists)

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
            labels: List<PeakAnnotation>,
            id: String,
            parameters: List<Pair<Double, Int>>,
            cancellableState: CancellableState = CancellableState.current()
    ): Pair<List<LabelErrors>, Int> {
        MultitaskProgress.addTask(id, parameters.size.toLong())
        val tuneResults = tune(results, results.fitInfo.genomeQuery(), labels, id, parameters, cancellableState)
        MultitaskProgress.finishTask(id)
        return tuneResults
    }


    fun tune(
            results: SpanFitResults,
            genomeQuery: GenomeQuery,
            labels: List<PeakAnnotation>,
            id: String,
            parameters: List<Pair<Double, Int>>,
            cancellableState: CancellableState = CancellableState.current()
    ): Pair<List<LabelErrors>, Int> {
        val labeledGenomeQuery = GenomeQuery(
            genomeQuery.genome,
            *labels.map { it.location.chromosome.name }.distinct().toTypedArray()
        )
        // Parallelism is OK here:
        // 1. getPeaks creates BitterSet per each parameters combination of size
        // ~ 3*10^9 / 200bp / 8 / 1024 / 1024 ~ 2MB for human genome
        // 2. List.parallelStream()....collect(Collectors.toList()) guarantees the same order as in original list.
        // Order is important!
        val labelErrorsGrid = Array<LabelErrors?>(parameters.size) { null }
        tuningExecutor.awaitAll(
                parameters.mapIndexed { index, (fdr, gap) ->
                    Callable {
                        cancellableState.checkCanceled()
                        val peaksOnLabeledGenomeQuery = results.getPeaks(labeledGenomeQuery, fdr, gap)
                        labelErrorsGrid[index] = computeErrors(labels,
                                LocationsMergingList.create(labeledGenomeQuery,
                                        peaksOnLabeledGenomeQuery.map { it.location }.iterator()))
                        MultitaskProgress.reportTask(id)
                    }
                }

        )
        val minTotalError = labelErrorsGrid.map { it!!.error() }.min()!!
        // Parameters should return desired order for each tool
        return labelErrorsGrid.map { it!! } to
                parameters.indices.first { labelErrorsGrid[it]!!.error() == minTotalError }
    }

}

object SpanReplicated : ReplicatedTool2Tune<Pair<Double, Int>>() {

    override val id = "span_replicated"
    override val suffix = Span.suffix
    const val DEFAULT_BIN = Span.DEFAULT_BIN

    private const val DEFAULT_FDR = 1E-60
    private val FDRS = doubleArrayOf(1E-10, 1E-20, 1E-30, 1E-40, 1E-50, DEFAULT_FDR, 1E-70, 1E-80, 1E-90, 1E-100)
    private val GAPS = Span.GAPS

    override val transform = Span.transform

    override val parameters =
            FDRS.sorted().flatMap { fdr ->
                GAPS.sorted().map { gap -> fdr to gap }
            }

    override fun defaultParams(uli: Boolean) = DEFAULT_FDR to Span.DEFAULT_GAP

    fun fileName(target: String, parameter: Pair<Double, Int>) =
            "${target}_${DEFAULT_BIN}_${parameter.first}_${parameter.second}_peaks.bed"

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
            PeakCallerTuning.LOG.info("$target $id Missing or extra peak files detected, launching tuning.")
            existingPeaks.forEach {
                PeakCallerTuning.LOG.info("Removing obsolete file $it")
                it.deleteIfExists()
            }
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

        val replicatedPeakCallingExperiment = SpanPeakCallingExperiment.getExperiment(
            configuration.genomeQuery,
            labelledTracks.map(LabelledTrack::trackPath).map { it to inputPath },
            DEFAULT_BIN, AutoFragment
        )

        if (saveAllPeaks) {
            PeakCallerTuning.LOG.info("Saving all $id peaks $target")
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

        PeakCallerTuning.LOG.info("Tuning $id peaks $target")
        val labelsPath = labelledTracks.first().labelPath
        val labels = PeakAnnotation.loadLabels(labelsPath, configuration.genomeQuery.genome)
        val (labelErrorsGrid, index) = Span.tune(replicatedPeakCallingExperiment.results, labels, target, parameters)

        PeakCallerTuning.LOG.info("Saving $target optimal $id peaks to $folder")
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