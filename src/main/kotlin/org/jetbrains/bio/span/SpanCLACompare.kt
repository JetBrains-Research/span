package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiment.configurePaths
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.span.SpanCLA.LOG
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanDifferentialPeakCallingExperiment
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.util.*
import org.slf4j.event.Level
import java.nio.file.Path

object SpanCLACompare {

    internal fun compare(params: Array<String>) {
        with(SpanCLA.getOptionParser()) {
            acceptsAll(
                listOf("t1", "treatment1"),
                """
                    ChIP-seq treatment file 1. bam, bed or .bed.gz file;
                    If multiple files are given, treated as replicates.
                    """.trimIndent()
            )
                .withRequiredArg().required()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("c1", "control1"),
                """
                    Control file 1. bam, bed or .bed.gz file;
                    Single control file or separate file per each
                    treatment file required.
                    """.trimIndent()
            )
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(
                listOf("t2", "treatment2"),
                """
                    ChIP-seq treatment file 2. bam, bed or .bed.gz file;
                    If multiple files are given, treated as replicates.
                    """.trimIndent()
            )
                .withRequiredArg().required()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("c2", "control2"),
                """
                    Control file 2. bam, bed or .bed.gz file;
                    Single control file or separate file per each
                    treatment file required.
                    """.trimIndent()
            )
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            accepts("clip", "Clip peaks to improve density")

            parse(params) { options ->
                if ("quiet" in options) {
                    Logs.quiet()
                } else {
                    Logs.addConsoleAppender(if ("debug" in options) Level.DEBUG else Level.INFO)
                }
                LOG.info("SPAN ${SpanCLA.version()}")
                LOG.info("COMMAND: compare ${params.joinToString(" ")}")

                SpanCLA.checkMemory()

                val peaksPath = options.valueOf("peaks") as Path?
                val labelsPath = options.valueOf("labels") as Path?
                val experimentId = if (peaksPath != null) {
                    peaksPath.stemGz
                } else {
                    SpanCLAAnalyze.generateExperimentId(
                        SpanCLAAnalyze.prepareAndCheckTreatmentControlPaths(options),
                        SpanCLA.getBin(options),
                        SpanCLA.getFragment(options),
                        SpanCLA.getUnique(options),
                        labelsPath
                    )
                }

                val workingDir = options.valueOf("workdir") as Path
                val logPath = options.valueOf("log") as Path?
                val chromSizesPath = options.valueOf("chrom.sizes") as Path?

                // Configure working directories
                LOG.info("WORKING DIR: $workingDir")
                if (!SpanCLA.ignoreConfigurePaths) {
                    configurePaths(workingDir, chromSizesPath = chromSizesPath, logPath = logPath)
                }
                // Configure logging to file
                val actualLogPath = logPath ?: (org.jetbrains.bio.experiment.Configuration.logsPath / "${experimentId}.log")
                Logs.addLoggingToFile(actualLogPath)
                LOG.info("LOG: $actualLogPath")


                // Call now to preserve correct params logging
                val lazyDifferentialPeakCallingResults = differentialPeakCallingResults(options)

                val fdr = options.valueOf("fdr") as Double
                val gap = options.valueOf("gap") as Int
                require(gap >= 0) { "Negative gap: $gap" }
                require(0 < fdr && fdr <= 1) { "Illegal fdr: $fdr, expected range: (0, 1)" }
                LOG.info("FDR: $fdr")
                LOG.info("GAP: $gap")

                if (peaksPath != null) {
                    LOG.info("PEAKS: $peaksPath")
                } else {
                    LOG.info("NO peaks path given, process model fitting only.")
                    LOG.info("LABELS, FDR, GAP options are ignored.")
                }

                val threads = options.valueOf("threads") as Int?
                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")

                val clip = "clip" in options
                LOG.info("CLIP: $clip")

                val differentialPeakCallingResults = lazyDifferentialPeakCallingResults.value
                val genomeQuery = differentialPeakCallingResults.fitInfo.genomeQuery()
                if (peaksPath != null) {
                    val peaks = ModelToPeaks.computeChromosomePeaks(
                        differentialPeakCallingResults,
                        genomeQuery,
                        fdr,
                        gap,
                        clip
                    )
                    Peak.savePeaks(peaks, peaksPath, "diff_${experimentId}_${fdr}_$gap")
                    LOG.info("Saved result to $peaksPath")
                }
                val keepCacheFiles = "keep-cache" in options
                if (!keepCacheFiles) {
                    LOG.debug("Clean coverage caches")
                    differentialPeakCallingResults.fitInfo.cleanCaches()
                }

            }
        }
    }


    /**
     * Retrieves the paths (treatment1, optional control1), (treatment2, optional control2)
     * Checks that they are consistent.
     */
    private fun getComparePaths(
        options: OptionSet,
        log: Boolean = false
    ): Pair<List<SpanDataPaths>, List<SpanDataPaths>> {
        val treatmentPaths1 = options.valuesOf("treatment1") as List<Path>
        val treatmentPaths2 = options.valuesOf("treatment2") as List<Path>
        val controlPaths1 = options.valuesOf("control1") as List<Path>
        val controlPaths2 = options.valuesOf("control2") as List<Path>

        val paths1 = SpanCLA.matchTreatmentsAndControls(treatmentPaths1, controlPaths1)
        check(paths1 != null) { "No treatment files provided for set 1, use -t1 option." }
        val paths2 = SpanCLA.matchTreatmentsAndControls(treatmentPaths2, controlPaths2)
        check(paths2 != null) { "No treatment files provided for set 2, use -t2 option." }

        if (log) {
            LOG.info("TREATMENT1: ${treatmentPaths1.joinToString(", ", transform = Path::toString)}")
            if (controlPaths1.isNotEmpty()) {
                LOG.info("CONTROL1: ${controlPaths1.joinToString(", ", transform = Path::toString)}")
            } else {
                LOG.info("CONTROL1: none")
            }
            LOG.info("TREATMENT2: ${treatmentPaths2.joinToString(", ", transform = Path::toString)}")
            if (controlPaths2.isNotEmpty()) {
                LOG.info("CONTROL2: ${controlPaths2.joinToString(", ", transform = Path::toString)}")
            } else {
                LOG.info("CONTROL2: none")
            }
        }
        return paths1 to paths2
    }

    /**
     * Configure logging and get [SpanFitResults] in a most concise and effective way.
     * Parses and logs most of the command line arguments.
     */
    private fun differentialPeakCallingResults(options: OptionSet): Lazy<SpanFitResults> {
        val workingDir = options.valueOf("workdir") as Path
        LOG.info("WORKING DIR: $workingDir")
        val chromSizesPath = options.valueOf("chrom.sizes") as Path?
        val genomeQuery = GenomeQuery(Genome[chromSizesPath!!])
        val (data1, data2) = getComparePaths(options, log = true)
        LOG.info("CHROM.SIZES: $chromSizesPath")
        val bin = SpanCLA.getBin(options, log = true)
        val fragment = SpanCLA.getFragment(options, log = true)
        val unique = SpanCLA.getUnique(options, log = true)
        val maxIterations = SpanCLA.getMaxIter(options, log = true)
        val threshold = SpanCLA.getThreshold(options, log = true)
        return lazy {
            val experiment = SpanDifferentialPeakCallingExperiment.getExperiment(
                genomeQuery, data1, data2, bin, fragment, unique,
                threshold, maxIterations
            )
            experiment.results
        }
    }

}