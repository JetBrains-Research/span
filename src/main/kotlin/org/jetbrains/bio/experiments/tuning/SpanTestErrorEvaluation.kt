package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.dataframe.DataFrameSpec
import org.jetbrains.bio.dataset.ChipSeqTarget
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.dataset.DataType
import org.jetbrains.bio.dataset.toDataType
import org.jetbrains.bio.experiment.Experiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.span.getPeaks
import org.jetbrains.bio.statistics.distribution.Sampling
import org.jetbrains.bio.util.div
import org.jetbrains.bio.util.toPath
import java.util.concurrent.ConcurrentHashMap

class SpanTestErrorEvaluation(
        val dataConfig: DataConfig, val ks: IntArray
) : Experiment("span-test-error") {

    private fun processModification(target: String): DataFrame {
        val builder = DataFrameSpec(synchronized = true)
                .strings("target", "cell")
                .ints("k", "batch")
                .doubles("train_error", "test_error", "fdr")
                .ints("gap")
                .builder()
        val tracks = dataConfig.extractLabelledTracks(target)
        tracks.forEach {
            println(it.name)
            val results = processTrack(target, it)
            ks.forEach { k ->
                results[k]?.forEachIndexed { batch, result ->
                    builder.add(target, it.name, k, batch, result.trainingError, result.testError, result.fdr, result.gap)
                }
            }
        }
        return builder.build()
    }

    internal data class ProcessingResults(val trainingError: Double, val testError: Double, val fdr: Double, val gap: Int)

    private fun processTrack(target: String, track: LabelledTrack): Map<Int, List<ProcessingResults>> {
        val name = track.name
        val labelPath = track.labelPath
        val labels = PeakAnnotation.loadLabels(labelPath, dataConfig.genomeQuery.genome)
        val trackPath = track.trackPath
        val inputPath = dataConfig.tracksMap.entries
                .filter { it.key.dataType.toDataType() == DataType.CHIP_SEQ && ChipSeqTarget.isInput(it.key.dataType) }
                .flatMap { it.value }.map { it.second.path }.first()
        val results = SpanPeakCallingExperiment.getExperiment(
            dataConfig.genomeQuery,
            listOf(Triple(trackPath, inputPath, inputPath)),
            Span.DEFAULT_BIN, AutoFragment
        ).results
        val res = ConcurrentHashMap<Int, List<ProcessingResults>>()
        ks.toList().parallelStream().forEach { k ->
            println("k: $k")
            val split = trainTestSplit(labels, k)
            val list = split.take(10).mapIndexed { batch, trainTest ->
                val (errors, optimal) = Span.tune(
                    results, trainTest.train, "$target-$name-$k-$batch", Span.parameters
                )
                val trainError = errors[optimal].error()
                val (fdr, gap) = Span.parameters[optimal]
                val peaks = results.getPeaks(dataConfig.genomeQuery, fdr, gap)
                val testError = computeErrors(
                    trainTest.test,
                    LocationsMergingList.create(dataConfig.genomeQuery, peaks.map { it.location })
                ).error()
                ProcessingResults(trainError, testError, fdr, gap)
            }
            res[k] = list
        }
        return res
    }

    override fun doCalculations() {
        dataConfig.dataTypes().filter {
            it.toDataType() == DataType.CHIP_SEQ && !ChipSeqTarget.isInput(it)
        }.forEach {
            processModification(it).save(experimentPath / "$it.tsv")
        }
    }

    data class TrainTest(val train: List<PeakAnnotation>, val test: List<PeakAnnotation>)

    companion object {

        fun trainTestSplit(labels: List<PeakAnnotation>, k: Int): List<TrainTest> {
            val numberOfBatches = labels.size / k
            if (numberOfBatches == 0) {
                return emptyList()
            }
            val numberOfLabels = numberOfBatches * k
            val availableLabels = Sampling.sampleCombination(labels.size, numberOfLabels).map { labels[it] }
            var remainingLabels = availableLabels
            return (0 until numberOfBatches).map {
                val trainingIndices = Sampling.sampleCombination(remainingLabels.size, k)
                val training = trainingIndices.map { remainingLabels[it] }
                val test = availableLabels.filter { it !in training }
                remainingLabels = remainingLabels.filter { it !in training }
                TrainTest(training, test)
            }
        }

        @JvmStatic
        fun main(args: Array<String>) {
            val dataConfig = DataConfig.load(args[0].toPath())
            val ks = args[1].split(",").map { Integer.parseInt(it) }.toIntArray()
            SpanTestErrorEvaluation(dataConfig, ks).run()
        }
    }
}
