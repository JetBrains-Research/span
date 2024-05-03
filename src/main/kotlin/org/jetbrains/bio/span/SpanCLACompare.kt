package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.experiment.configurePaths
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.span.SpanCLA.LOG
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BACKGROUND_SENSITIVITY
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
                    ChIP-seq treatment file 1. bam, bed or .bed.gz file.
                    If multiple files are given, treated as replicates
                    """.trimIndent()
            )
                .withRequiredArg().required()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("c1", "control1"),
                """
                    Control file 1. bam, bed or .bed.gz file.
                    Single control file or separate file per each treatment file required
                    """.trimIndent()
            )
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(
                listOf("t2", "treatment2"),
                """
                    ChIP-seq treatment file 2. bam, bed or .bed.gz file.
                    If multiple files are given, treated as replicates
                    """.trimIndent()
            )
                .withRequiredArg().required()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("c2", "control2"),
                """
                    Control file 2. bam, bed or .bed.gz file.
                    Single control file or separate file per each treatment file required
                    """.trimIndent()
            )
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            accepts("clip", "Clip peaks to improve peaks density using local signal coverage.\n" +
                    "Recommended for TFs, narrow histone marks and ATAC-seq")

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
                val keepCacheFiles = "keep-cache" in options
                checkOrFail(peaksPath != null || keepCacheFiles) {
                    "At least one of the parameters is required: --peaks or --keep-cache."
                }

                val modelId = peaksPath?.stemGz ?:
                SpanAnalyzeFitInformation.generateId(
                    SpanCLAAnalyze.prepareAndCheckTreatmentControlPaths(options),
                    SpanCLA.getFragment(options),
                    SpanCLA.getBin(options),
                    SpanCLA.getUnique(options),
                )

                val fdr = options.valueOf("fdr") as Double
                require(0 < fdr && fdr < 1) { "Illegal fdr: $fdr, expected range: (0, 1)" }
                val bgSensitivity = if (options.has("bg-sensitivity"))
                    options.valueOf("bg-sensitivity") as Double
                else
                    SPAN_DEFAULT_BACKGROUND_SENSITIVITY
                require(bgSensitivity.isNaN() || 0 < bgSensitivity && bgSensitivity <= 1) {
                    "Illegal background sensitivity: $bgSensitivity, expected range: (0, 1]"
                }
                val clip = if (options.has("clip")) options.valueOf("clip") as Double else SPAN_DEFAULT_CLIP
                require(0 <= clip && clip < 1) { "Illegal clip: $fdr, expected range: [0, 1)" }

                val workingDir = options.valueOf("workdir") as Path
                val id = peaksPath?.stemGz ?: reduceIds(listOf(modelId, fdr.toString(), bgSensitivity.toString(), clip.toString()))
                var logPath = options.valueOf("log") as Path?
                val chromSizesPath = options.valueOf("chrom.sizes") as Path?

                // Configure working directories
                LOG.info("WORKING DIR: $workingDir")
                if (!SpanCLA.ignoreConfigurePaths) {
                    configurePaths(workingDir, genomesPath = chromSizesPath?.parent, logsPath = logPath?.parent)
                }
                // Configure logging to file
                if (logPath == null) {
                    logPath = Configuration.logsPath / "$id.log"
                }
                Logs.addLoggingToFile(logPath)
                LOG.info("LOG: $logPath")

                // Call now to preserve correct params logging
                val lazyDifferentialPeakCallingResults = differentialPeakCallingResults(options)

                LOG.info("FDR: $fdr")
                LOG.info("BACKGROUND SENSITIVITY: $bgSensitivity")
                LOG.info("CLIP: $clip")

                if (peaksPath != null) {
                    LOG.info("PEAKS: $peaksPath")
                } else {
                    LOG.info("NO peaks path given, process model fitting only.")
                    LOG.info("Labels, fdr, background sensitivity, clip options are ignored.")
                }

                val threads = options.valueOf("threads") as Int? ?: Runtime.getRuntime().availableProcessors()
                check(threads > 0) {
                    "Negative threads value: $threads"
                }
                check(threads <= Runtime.getRuntime().availableProcessors()) {
                    "Too big threads value $threads > ${Runtime.getRuntime().availableProcessors()}"
                }
                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")


                val differentialPeakCallingResults = lazyDifferentialPeakCallingResults.value
                val genomeQuery = differentialPeakCallingResults.fitInfo.genomeQuery()
                if (peaksPath != null) {
                    val peaks = ModelToPeaks.getPeaks(
                        differentialPeakCallingResults,
                        genomeQuery,
                        fdr,
                        bgSensitivity,
                        clip
                    )
                    Peak.savePeaks(peaks, peaksPath, "diff_${id}.peak")
                    LOG.info("Saved result to $peaksPath")
                }
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