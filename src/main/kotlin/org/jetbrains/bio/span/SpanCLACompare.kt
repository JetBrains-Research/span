package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.experiment.configurePaths
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.span.SpanCLA.LOG
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_HARD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_SPEED
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_LIGHT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION
import org.jetbrains.bio.span.fit.SpanConstants.printSpanConstants
import org.jetbrains.bio.span.peaks.SpanModelToPeaks
import org.jetbrains.bio.span.peaks.MultipleTesting
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
                val sensitivity = if (options.has("sensitivity")) options.valueOf("sensitivity") as Double else null
                val gap = if (options.has("gap")) options.valueOf("gap") as Int else null
                require(gap == null || gap >= 0) { "Illegal gap: $gap, expected >= 0" }

                val workingDir = options.valueOf("workdir") as Path
                val id = peaksPath?.stemGz ?:
                reduceIds(listOfNotNull(modelId, fdr.toString(), sensitivity?.toString(), gap?.toString()))
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

                val clip = if (options.has("clip")) options.valueOf("clip") as Double else SPAN_DEFAULT_CLIP_MAX_SIGNAL

                if (peaksPath != null) {
                    LOG.info("FDR: $fdr")
                    LOG.info("SENSITIVITY: $sensitivity")
                    LOG.info("GAP: $gap")
                    LOG.info("CLIP: $clip")
                    LOG.info("PEAKS: $peaksPath")
                } else {
                    LOG.info("NO peaks path given, process model fitting only.")
                    LOG.info("Labels, fdr, sensitivity, gap, clip options are ignored.")
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

                // Print all the constants, which are not configured using command line
                if (LOG.isDebugEnabled) {
                    printSpanConstants()
                }

                val differentialPeakCallingResults = lazyDifferentialPeakCallingResults.value
                val genomeQuery = differentialPeakCallingResults.fitInfo.genomeQuery()

                val multipleTesting = if ("multiple" in params)
                    MultipleTesting.valueOf(options.valueOf("multiple") as String)
                else
                    SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION
                LOG.info("MULTIPLE TEST CORRECTION: ${multipleTesting.description}")

                val blackListPath = options.valueOf("blacklist") as Path?

                if (peaksPath != null) {
                    val peaks = SpanModelToPeaks.getPeaks(
                        differentialPeakCallingResults,
                        genomeQuery,
                        fdr, multipleTesting,
                        sensitivity, gap, false,
                        SPAN_DEFAULT_FRAGMENTATION_LIGHT,
                        SPAN_DEFAULT_FRAGMENTATION_HARD,
                        SPAN_DEFAULT_FRAGMENTATION_SPEED,
                        clip = clip,
                        blackListPath = blackListPath,
                        name = peaksPath.fileName.stem
                    )
                    LOG.info("Format chromosome, start, end, name, score, strand, foldchange, -log(p), -log(q)")
                    Peak.savePeaks(peaks.toList(), peaksPath, "diff_${id}.peak")
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
        val fitThreshold = SpanCLA.getFitThreshold(options, log = true)
        val fitMaxIterations = SpanCLA.getFitMaxIteration(options, log = true)
        return lazy {
            val experiment = SpanDifferentialPeakCallingExperiment.getExperiment(
                genomeQuery, data1, data2, bin, fragment, unique,
                fitThreshold, fitMaxIterations
            )
            experiment.results
        }
    }

}