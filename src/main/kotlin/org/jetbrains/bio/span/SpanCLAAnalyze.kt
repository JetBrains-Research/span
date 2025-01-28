package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.experiment.configurePaths
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.SpanCLA.LOG
import org.jetbrains.bio.span.SpanCLA.checkGenomeInFitInformation
import org.jetbrains.bio.span.SpanResultsAnalysis.doDeepAnalysis
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_COMPENSATION_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_LOW_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_HMM_ESTIMATE_SNR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_MAX_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.printSpanConstants
import org.jetbrains.bio.span.fit.experimental.*
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.peaks.MultipleTesting
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.peaks.SpanPeaksResult
import org.jetbrains.bio.span.semisupervised.LocationLabel
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised
import org.jetbrains.bio.span.semisupervised.TuningResults
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.util.*
import org.slf4j.event.Level
import java.nio.file.Path


object SpanCLAAnalyze {

    internal fun analyze(params: Array<String>) {
        with(SpanCLA.getOptionParser()) {
            acceptsAll(
                listOf("t", "treatment"),
                """
                    ChIP-seq treatment file. bam, bed or .bed.gz file.
                    If multiple files are given, treated as replicates
                    """.trimIndent()
            )
                .requiredUnless("model")
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("c", "control"),
                """
                    Control file. bam, bed or bed.gz file.
                    Single control file or separate file per each treatment file required
                    """.trimIndent()
            )
                .availableIf("treatment")
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.noCheck())

            accepts("bigwig", "Create beta-control corrected counts per million normalized track")
                .availableIf("peaks")

            accepts("labels", "Labels BED file. Used in semi-supervised peak calling")
                .availableIf("peaks")
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.exists())
            accepts(
                "model-type",
                "Model type. Experimental.\n" +
                        SpanModelType.values().joinToString("\n") { "${it.id}\t${it.description}" }
            )
                .withRequiredArg()
                .defaultsTo(SpanModelType.NB2Z_HMM.id)

            accepts(
                "ext", "Save extended states information to model file.\n" +
                        "Required for model visualization in JBR Genome Browser"
            )

            accepts(
                "deep-analysis",
                "Deep analysis of model including analysis of coverage / candidates / peaks"
            )

            accepts("hmm-snr",
                "Fraction of coverage to estimate and guard signal to noise ratio, 0 to disable")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_HMM_ESTIMATE_SNR)

            accepts("hmm-low",
                "Minimal low state mean threshold, guards against too broad peaks, 0 to disable")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_HMM_LOW_THRESHOLD)

            accepts("fragmentation",
                "Fragmentation threshold to enable compensation")
                .availableUnless("gap")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_FRAGMENTATION_MAX_THRESHOLD)


            accepts("gap-fragmented",
                "Additional compensation gap for tracks with high fragmentation")
                .availableUnless("gap")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(SPAN_DEFAULT_FRAGMENTATION_COMPENSATION_GAP)

            // Additional parameter for *_REGRESSION_MIXTURE models
            accepts(
                "mapability",
                "Mapability bigWig file. Requires --model-type " +
                        "${SpanModelType.POISSON_REGRESSION_MIXTURE.id} or " +
                        "--model-type ${SpanModelType.NEGBIN_REGRESSION_MIXTURE.id}"

            )
                .availableIf("treatment")
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->
                if ("quiet" in options) {
                    Logs.quiet()
                } else {
                    Logs.addConsoleAppender(if ("debug" in options) Level.DEBUG else Level.INFO)
                }
                LOG.info("SPAN ${SpanCLA.version()}")
                LOG.info(
                    "COMMAND: analyze ${params.joinToString(" ")}"
                )

                SpanCLA.checkMemory()
                val peaksPath = options.valueOf("peaks") as Path?
                val modelPath = options.valueOf("model") as Path?
                val deepAnalysis = "deep-analysis" in options
                val blackListPath = options.valueOf("blacklist") as Path?
                val keepCacheFiles = "keep-cache" in options
                checkOrFail(peaksPath != null || modelPath != null || keepCacheFiles) {
                    "At least one of the parameters is required: --peaks, --model or --keep-cache."
                }

                val labelsPath = options.valueOf("labels") as Path?
                val modelId = peaksPath?.stemGz ?: modelPath?.stem ?: SpanAnalyzeFitInformation.generateId(
                    prepareAndCheckTreatmentControlPaths(options),
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
                val id = peaksPath?.stemGz ?: reduceIds(
                    listOfNotNull(
                        modelId,
                        fdr.toString(),
                        sensitivity?.toString(),
                        gap?.toString()
                    )
                )
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

                // Call now to preserve params logging order
                val lazySpanResults = logParametersAndPrepareLazySpanResults(options)

                val bigWig = options.contains("bigwig")
                LOG.info("BIGWIG: $bigWig")

                val multipleTesting = if ("multiple" in params)
                    MultipleTesting.valueOf(options.valueOf("multiple") as String)
                else
                    SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION
                LOG.info("MULTIPLE TEST CORRECTION: ${multipleTesting.description}")

                val clip = options.valueOf("clip") as Double
                val fragmentationThreshold  = when {
                    gap != null -> 0.0
                    options.has("fragmentation") -> options.valueOf("fragmentation") as Double
                    else -> SPAN_DEFAULT_FRAGMENTATION_MAX_THRESHOLD
                }

                val gapFragmentationCompensation = when {
                    gap != null -> 0
                    options.has("gap-fragmentation") -> options.valueOf("gap-fragmentation") as Int
                    else -> SPAN_DEFAULT_FRAGMENTATION_COMPENSATION_GAP
                }

                if (peaksPath != null) {
                    if (labelsPath != null) {
                        LOG.info("LABELS: $labelsPath")
                        LOG.info("Fdr, sensitivity, gap, and other options are ignored.")
                    } else {
                        LOG.info("FDR: $fdr")
                        if (sensitivity != null) {
                            LOG.info("SENSITIVITY: $sensitivity")
                        }
                        if (gap != null) {
                            LOG.info("GAP: $gap")
                        }
                    }
                    if (gap == null) {
                        LOG.info("FRAGMENTATION THRESHOLD: $fragmentationThreshold")
                        LOG.info("GAP FRAGMENTATION: $gapFragmentationCompensation")
                    }
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

                // Finally get SPAN results
                val (actualModelPath, spanResults) = lazySpanResults.value
                val fitInfo = spanResults.fitInfo
                check(fitInfo is AbstractSpanAnalyzeFitInformation) {
                    "Expected ${SpanAnalyzeFitInformation::class.java.simpleName}, got ${fitInfo::class.java.name}"
                }
                val aboutModel = spanResults.modelInformation(actualModelPath)
                LOG.info(aboutModel.joinToString("\n") { (k, v) ->
                    "${k.name}: ${k.render(v)}"
                })
                LOG.debug(spanResults.model.toString())

                val genomeQuery = fitInfo.genomeQuery()
                val fragment = fitInfo.fragment
                val bin = fitInfo.binSize

                if (bigWig) {
                    if (fitInfo !is SpanAnalyzeFitInformation) {
                        LOG.warn("Bigwig coverage is possible only for analyze command")
                    } else {
                        val bigWigPath = (peaksPath!!.toString() + ".bw").toPath()
                        SpanBigWigWriter.write(spanResults, genomeQuery, bigWigPath, blackListPath)
                    }
                }

                if (deepAnalysis) {
                    if (fitInfo !is SpanAnalyzeFitInformation) {
                        LOG.warn("Deep analysis is possible only for analyze command")
                    }
                }
                if (peaksPath != null) {
                    val peaks = if (labelsPath == null)
                        ModelToPeaks.getPeaks(
                            spanResults, genomeQuery, fdr, multipleTesting,
                            sensitivity, gap,
                            fragmentationThreshold, gapFragmentationCompensation,
                            clip = clip,
                            blackListPath = blackListPath
                        )
                    else
                        tune(spanResults, genomeQuery, labelsPath, peaksPath, blackListPath)
                    val peaksList = processBlackList(genomeQuery, peaks, blackListPath)

                    LOG.info("Format chromosome, start, end, name, score, strand, signal, -log(p), -log(q)")
                    Peak.savePeaks(
                        peaksList, peaksPath,
                        "peak${if (fragment is FixedFragment) "_$fragment" else ""}_" +
                                "${bin}_${fdr}_${peaks.sensitivity}_${peaks.gap}"
                    )
                    LOG.info("Peaks saved to $peaksPath")
                    val aboutPeaks = PeaksInfo.compute(
                        genomeQuery,
                        peaksList.map { it.location }.stream(),
                        peaksPath.toUri(),
                        fitInfo.paths.map { it.treatment }
                    )
                    LOG.info(aboutPeaks.joinToString("\n") { (k, v) ->
                        "${k.name}: ${k.render(v)}"
                    })

                    if (deepAnalysis && fitInfo is SpanAnalyzeFitInformation) {
                        fitInfo.prepareData()
                        doDeepAnalysis(
                            actualModelPath,
                            spanResults,
                            fitInfo,
                            genomeQuery,
                            fdr,
                            sensitivity, gap,
                            fragmentationThreshold,
                            gapFragmentationCompensation,
                            blackListPath,
                            peaksList,
                            peaksPath
                        )
                    }
                }
                if (modelPath == null && !keepCacheFiles) {
                    LOG.debug("Clean coverage caches")
                    fitInfo.cleanCaches()
                }
            }
        }
    }

    fun processBlackList(
        genomeQuery: GenomeQuery,
        peaks: SpanPeaksResult,
        blacklistPath: Path?
    ): List<Peak> {
        var peaksList = peaks.toList()
        if (blacklistPath != null) {
            LOG.info("Filter out blacklisted regions")
            val blackList = LocationsMergingList.load(genomeQuery, blacklistPath)
            peaksList = peaksList.filter { !blackList.intersects(it.location) }
        }
        return peaksList
    }


    private fun tune(
        spanResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        labelsPath: Path,
        peaksPath: Path,
        blackListPath: Path? = null,
    ): SpanPeaksResult {
        val results = TuningResults()
        LOG.info("Loading labels $labelsPath...")
        val labels = LocationLabel.loadLabels(labelsPath, genomeQuery.genome)
        LOG.info("Tuning model on the loaded labels...")
        val (labelErrorsGrid, optimalIndex) = SpanSemiSupervised.tuneParameters(
            spanResults,
            genomeQuery,
            labels,
            "",
            SpanSemiSupervised.PARAMETERS,
            CancellableState.current()
        )
        LOG.info("Tuning model on the loaded labels complete.")
        val (optimalFDR, optimalSensitivity, optimalGap) = SpanSemiSupervised.PARAMETERS[optimalIndex]
        LOG.info("Optimal settings: FDR=$optimalFDR, SENSITIVITY=$optimalSensitivity, GAP=$optimalGap")
        labelErrorsGrid.forEachIndexed { i, error ->
            val (fdrTuning, sensitivityTuning, gapTuning) = SpanSemiSupervised.PARAMETERS[i]
            results.addRecord(
                "result",
                "${fdrTuning}_${sensitivityTuning}_${gapTuning}",
                error,
                i == optimalIndex
            )
        }
        results.saveTuningErrors(peaksPath.parent / "${peaksPath.fileName.stem}_errors.csv")
        results.saveOptimalResults(
            peaksPath.parent
                    / "${peaksPath.fileName.stem}_parameters.csv"
        )
        return ModelToPeaks.getPeaks(
            spanResults, genomeQuery, optimalFDR, SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION,
            optimalSensitivity, optimalGap,
            SPAN_DEFAULT_FRAGMENTATION_MAX_THRESHOLD,
            SPAN_DEFAULT_FRAGMENTATION_COMPENSATION_GAP,
            SPAN_DEFAULT_CLIP_MAX_SIGNAL,
            blackListPath = blackListPath
        )
    }

    /**
     * Retrieves the paths (treatment, optional control)
     * either from command-line options or from the stored fit information.
     *
     * If both are available, checks that they are consistent.
     */
    internal fun prepareAndCheckTreatmentControlPaths(
        options: OptionSet, fitInformation: AbstractSpanAnalyzeFitInformation? = null, log: Boolean = false
    ): List<SpanDataPaths> {
        val commandLineTreatmentPaths = options.valuesOf("treatment") as List<Path>
        val commandLineControlPaths = options.valuesOf("control") as List<Path>

        var paths = SpanCLA.matchTreatmentsAndControls(
            commandLineTreatmentPaths, commandLineControlPaths
        )

        val fitInfoPaths = fitInformation?.paths
        if (fitInfoPaths != null) {
            if (paths != null) {
                check(paths == fitInfoPaths) {
                    "Stored treatment-control pairs ${fitInfoPaths.joinToString()} differ from the ones inferred " +
                            "from the command line arguments: ${paths!!.joinToString()}"
                }
            } else {
                paths = fitInfoPaths
            }
        }

        checkNotNull(paths) {
            "No treatment files and no existing model file provided, exiting."
        }

        if (log) {
            LOG.info("TREATMENT: ${paths.map { it.treatment }.joinToString(", ", transform = Path::toString)}")
            paths.mapNotNull { it.control }.let {
                if (it.isNotEmpty()) {
                    LOG.info("CONTROL: ${it.joinToString(", ", transform = Path::toString)}")
                } else {
                    LOG.info("CONTROL: none")
                }
            }
        }
        return paths
    }

    /**
     * Log parameters and optionally computation of SPAN model, represented by [SpanFitResults].
     * Doesn't fit the model; this happens only if .value is invoked
     */
    private fun logParametersAndPrepareLazySpanResults(options: OptionSet): Lazy<Pair<Path, SpanFitResults>> {
        val modelPath = options.valueOf("model") as Path?
        if (modelPath != null) {
            LOG.info("MODEL: $modelPath")
        }
        if (modelPath != null && modelPath.exists && modelPath.size.isNotEmpty()) {
            LOG.debug(
                "Model file {} exists and is not empty, SPAN will use it to substitute the missing " +
                        "command line arguments and verify the provided ones.",
                modelPath
            )
            val results = SpanModelFitExperiment.loadResults(modelPath = modelPath)
            check(results.fitInfo is AbstractSpanAnalyzeFitInformation) {
                "Expected ${SpanAnalyzeFitInformation::class.java.simpleName}, got ${results.fitInfo::class.java.name}"
            }
            prepareAndCheckTreatmentControlPaths(options, results.fitInfo, log = true)
            val blacklistPath = options.valueOf("blacklist") as Path?
            LOG.info("BLACKLIST FILE: $blacklistPath")
            val workingDir = options.valueOf("workdir") as Path
            LOG.info("WORKING DIR: $workingDir")
            val chromSizesPath = options.valueOf("chrom.sizes") as Path?
            LOG.info("CHROM.SIZES: $chromSizesPath")
            val chromosomesToProcess = if (options.has("chromosomes"))
                options.valueOf("chromosomes").toString().split(',', ' ')
            else
                null
            LOG.info("CHROMOSOMES: ${chromosomesToProcess?.joinToString(", ")}")
            if (chromSizesPath != null) {
                checkGenomeInFitInformation(chromSizesPath, results.fitInfo)
            }
            SpanCLA.getBin(options, results.fitInfo, log = true)
            SpanCLA.getFragment(options, results.fitInfo, log = true)
            SpanCLA.getUnique(options, results.fitInfo, log = true)
            if (results.fitInfo is SpanRegrMixtureAnalyzeFitInformation) {
                getMapabilityPath(options, results.fitInfo, log = true)
            }
            val keepCacheFiles = "keep-cache" in options
            LOG.info("KEEP-CACHE: $keepCacheFiles")
            if (!keepCacheFiles) {
                LOG.warn("Keep cache files setting (false) is discarded, since fitted model already exists.")
            }
            return lazyOf(modelPath to results)        // Create fake lazy of already computed results
        } else {
            val paths = prepareAndCheckTreatmentControlPaths(options, log = true)
            val blacklistPath = options.valueOf("blacklist") as Path?
            LOG.info("BLACKLIST FILE: $blacklistPath")
            val workingDir = options.valueOf("workdir") as Path
            LOG.info("WORKING DIR: $workingDir")
            val chromSizesPath = options.valueOf("chrom.sizes") as Path?
            LOG.info("CHROM.SIZES: $chromSizesPath")
            val chromosomesToProcess = if (options.has("chromosomes"))
                options.valueOf("chromosomes").toString().split(',', ' ')
            else
                null
            LOG.info("CHROMOSOMES: ${chromosomesToProcess?.joinToString(", ")}")
            val bin = SpanCLA.getBin(options, log = true)
            val fragment = SpanCLA.getFragment(options, log = true)
            val unique = SpanCLA.getUnique(options, log = true)
            val fitThreshold = SpanCLA.getFitThreshold(options, log = true)
            val fitMaxIterations = SpanCLA.getFitMaxIteration(options, log = true)
            val mapabilityPath: Path?
            val modelType: SpanModelType = getModelType(options, modelPath)
            mapabilityPath = if (
                modelType == SpanModelType.POISSON_REGRESSION_MIXTURE ||
                modelType == SpanModelType.NEGBIN_REGRESSION_MIXTURE
            ) {
                getMapabilityPath(options, log = true)
            } else {
                null
            }
            val hmmEstimateSNR = options.valueOf("hmm-snr") as Double
            LOG.info("HMM ESTIMATE SNR: $hmmEstimateSNR")
            val hmmLow = options.valueOf("hmm-low") as Double
            LOG.info("HMM LOW STATE MIN: $hmmLow")
            val saveExtendedInfo = options.has("ext")
            LOG.info("EXTENDED MODEL INFO: $saveExtendedInfo")
            val keepCacheFiles = "keep-cache" in options
            LOG.info("KEEP-CACHE: $keepCacheFiles")
            return lazy {
                val genomeQuery = if (chromosomesToProcess != null)
                    GenomeQuery(Genome[chromSizesPath!!], *chromosomesToProcess.toTypedArray())
                else
                    GenomeQuery(Genome[chromSizesPath!!])
                val experiment = getExperimentByModelType(
                    modelType,
                    genomeQuery, paths, unique, fragment, bin,
                    hmmEstimateSNR, hmmLow,
                    fitThreshold, fitMaxIterations,
                    modelPath, saveExtendedInfo, keepCacheFiles,
                    mapabilityPath
                )
                experiment.modelPath to experiment.results
            }
        }
    }

    fun getExperimentByModelType(
        modelType: SpanModelType,
        genomeQuery: GenomeQuery,
        paths: List<SpanDataPaths>,
        unique: Boolean,
        fragment: Fragment,
        bin: Int,
        hmmEstimateSNR: Double,
        hmmLow: Double,
        fitThreshold: Double,
        fitMaxIterations: Int,
        modelPath: Path?,
        saveExtendedInfo: Boolean,
        keepCacheFiles: Boolean,
        mapabilityPath: Path?
    ): SpanModelFitExperiment<ClassificationModel, AbstractSpanAnalyzeFitInformation, out Enum<*>> {
        val experiment = when (modelType) {
            SpanModelType.NB2Z_HMM ->
                SpanPeakCallingExperiment.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, hmmEstimateSNR, hmmLow,
                    modelPath, fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB2Z_MIXTURE ->
                SpanPeakCallingExperimentNB2ZMixture.getExperiment(
                    genomeQuery, paths, fragment, bin, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB2_HMM ->
                SpanPeakCallingExperimentNB2HMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB3Z_HMM ->
                SpanPeakCallingExperimentNB3ZHMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB5Z_HMM ->
                SpanPeakCallingExperimentNB5ZHMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB3_HMM ->
                SpanPeakCallingExperimentNB3HMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.POISSON_REGRESSION_MIXTURE -> {
                SpanPeakCallingExperimentP2ZRegrMixture.getExperiment(
                    genomeQuery, paths, mapabilityPath, fragment, bin, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )
            }

            SpanModelType.NEGBIN_REGRESSION_MIXTURE -> {
                SpanPeakCallingExperimentNB2ZRegrMixture.getExperiment(
                    genomeQuery, paths, mapabilityPath, fragment, bin, unique, modelPath,
                    fitThreshold, fitMaxIterations, saveExtendedInfo, keepCacheFiles
                )
            }
        }
        return experiment
    }

    private fun getMapabilityPath(
        options: OptionSet, fitInfo: SpanRegrMixtureAnalyzeFitInformation? = null, log: Boolean = false
    ) = SpanCLA.getProperty(
        options.valueOf("mapability") as Path?,
        fitInfo?.mapabilityPath,
        null,
        "mapability", "MAPABILITY",
        log
    )

    private fun getModelType(
        options: OptionSet, modelPath: Path?
    ) = SpanCLA.getProperty(
        options.valueOf("model-type")?.let { SpanModelType.guessSpanModelById(it as String) },
        modelPath?.let { SpanModelType.guessSpanModelByExtension(it.extension) },
        SpanModelType.NB2Z_HMM, "model type", "MODEL TYPE", true
    )

}
