package org.jetbrains.bio.experiments.tuning

import joptsimple.OptionParser
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.Logs
import org.jetbrains.bio.datasets.ChipSeqTarget
import org.jetbrains.bio.datasets.DataType
import org.jetbrains.bio.datasets.toDataType
import org.jetbrains.bio.experiments.DataConfig
import org.jetbrains.bio.experiments.DataConfigExperiment
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.tools.Washu
import org.jetbrains.bio.util.PathConverter
import org.jetbrains.bio.util.contains
import org.jetbrains.bio.util.div
import java.nio.file.Path
import java.nio.file.Paths

/**
 * Created by Aleksei Dievskii on 24.10.2017.
 */

class PeakCallerTuning(configuration: DataConfig,
                       washuPath: Path = Washu.PATH,
                       val tools: List<Tool2Tune<*>>,
                       val useInput: Boolean = true,
                       private val computeTestError: Boolean = false) : DataConfigExperiment("benchmark", configuration) {

    val washu = Washu(washuPath)

    override fun doCalculations() {
        LOG.info("Processing peak caller tuning")
        val targets = configuration.dataTypes()
                .filter { it.toDataType() == DataType.CHIP_SEQ && !ChipSeqTarget.isInput(it) }
        targets.forEach { target ->
            LOG.info("Processing target $target")
            tools.forEach { t ->
                if (ChipSeqTarget.isWide(target) && t == MACS2) {
                    LOG.info("${MACS2.id} doesn't support wide modification $target")
                } else if (!ChipSeqTarget.isWide(target) && t == SICER) {
                    LOG.info("${SICER.id} doesn't support narrow modification $target")
                } else {
                    t.tune(configuration, experimentPath, target, useInput, true)
                    if (computeTestError) {
                        computeTestError(target, t)
                    }
                    generateConsensus(target, t)
                    if (t is ReplicatedTool2Tune) {
                        generateReplicatedLabelErrors(target, t)
                    } else {
                        generateLabelErrors(target, t)
                    }
                }
            }
        }

        LOG.info("Creating reports")
        val report = report()
        val uliReport = report()
        targets.forEach { target ->
            LOG.info("Processing target $target")
            tools.forEach { t ->
                if (ChipSeqTarget.isWide(target) && t == MACS2) {
                    LOG.info("${MACS2.id} doesn't support wide modification $target")
                } else if (!ChipSeqTarget.isWide(target) && t == SICER) {
                    LOG.info("${SICER.id} doesn't support narrow modification $target")
                } else {
                    computeFripAndReport(report, target, t,
                            t.folder(experimentPath, target, useInput), "tuned", washu)
                    computeFripAndReport(report, target, t,
                            t.defaultsFolder(experimentPath, target, useInput, false),
                            "default", washu)

                    computeFripAndReport(uliReport, target, t,
                            t.folder(experimentPath, target, useInput), "tuned", washu)
                    computeFripAndReport(uliReport, target, t,
                            t.defaultsFolder(experimentPath, target, useInput, true),
                            "default", washu)
                }
            }
        }
        val reportPath = experimentPath / "${configuration.id}_peaks_summary.tsv"
        report.build().save(reportPath)
        LOG.info("Defaults summary saved to $reportPath")
        val uliReportPath = experimentPath / "${configuration.id}_peaks_summary_uli.tsv"
        uliReport.build().save(uliReportPath)
        LOG.info("ULI defaults summary saved to $uliReportPath")
    }


    private fun computeTestError(target: String, tool: Tool2Tune<*>) {
        val tunedPeaks = tool.tunedPeaks(configuration, experimentPath, target, useInput)
        val labelledTracks = configuration.extractLabelledTracks(target).filter { it.testLabelPath != null }
        val resultTable = TuningResultTable()
        labelledTracks.forEach { track ->
            val peakPath = tunedPeaks[track.cellId to track.name]
            if (peakPath == null) {
                LOG.warn("Skipping non-existent peak file for $target $tool ${track.cellId} ${track.name}")
                return@forEach
            }
            resultTable.addRecord(track.name, "test",
                    computeErrors(PeakAnnotation.loadLabels(track.testLabelPath!!, configuration.genome),
                            LocationsMergingList.load(configuration.genomeQuery, peakPath)))
        }
        resultTable.print(tool.folder(experimentPath, target, useInput) / "${target}_${tool}_test_error.tsv")
    }


    companion object {
        internal val LOG = Logger.getLogger(PeakCallerTuning::class.java)

        @JvmStatic
        fun main(args: Array<String>) {
            Logs.addConsoleAppender(Level.INFO)
            val options = OptionParser().apply {
                nonOptions("config file")
                accepts("washu", "Washu scripts path")
                        .withRequiredArg().defaultsTo(Washu.PATH.toString())
                        .withValuesConvertedBy(PathConverter.exists())
                accepts("dir", "Working dir").withRequiredArg()
                        .defaultsTo(Paths.get(".").toAbsolutePath().toString()) // curr dir
                        .withValuesConvertedBy(PathConverter.exists())
                accepts("p", "Run SPAN")
                accepts("r", "Run replicated SPAN")
                accepts("s", "Run SICER")
                accepts("b", "Run MACS2 broad")
                accepts("n", "Run MACS2 narrow")
                accepts("i", "Don't use input (currently supported for SPAN only)")
                accepts("d", "Turn on debugging output")
                accepts("t", "Test result on a specified label file")
                acceptsAll(listOf("h", "?", "help"), "Show help").forHelp()
            }.parse(*args)

            if ("d" in options) {
                Logger.getRootLogger().level = Level.DEBUG
            }

            if (options.nonOptionArguments().isEmpty()) {
                LOG.error("No configuration given, nothing to run.")
                return
            }
            val config = options.nonOptionArguments()[0] as String
            val washuPath = if (options.has("washu"))
                options.valueOf("washu") as Path
            else
                Washu.PATH

            val tools = arrayListOf<Tool2Tune<*>>()
            if (options.has("p")) {
                tools.add(SPAN)
            }
            if (options.has("r")) {
                tools.add(SPAN_REPLICATED)
            }

            if (options.has("n")) {
                tools.add(MACS2)
            }
            if (options.has("b")) {
                tools.add(MACS2_BROAD)
            }
            if (options.has("s")) {
                tools.add(SICER)
            }

            val input = !options.has("i")
            val test = options.has("t")

            // Configuration
            Configuration.setExperimentWorkingDir(options.valueOf("dir") as Path)

            LOG.info("CONFIG:\t$config")
            LOG.info("WASHU PATH:\t$washuPath")
            LOG.info("WORKDIR:\t${Configuration.experimentsPath}")
            LOG.info("TOOLS:\t${tools.joinToString(",") { it.id }}")
            LOG.info("Use input: $input")
            LOG.info("Compute test error: $test")

            if (tools.isEmpty()) {
                LOG.info("No peak callers selected, nothing to run.")
            } else {
                PeakCallerTuning(loadDataConfig(config), washuPath,
                        tools = tools,
                        useInput = input,
                        computeTestError = test).run()
            }
        }
    }

}