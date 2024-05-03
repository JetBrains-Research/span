package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.experiment.configurePaths
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.SpanCLA.LOG
import org.jetbrains.bio.span.SpanCLA.checkGenomeInFitInformation
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP
import org.jetbrains.bio.span.fit.experimental.*
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.peaks.Peak
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
                val bgSensitivity = if (options.has("bg-sensitivity"))
                    options.valueOf("bg-sensitivity") as Double
                else
                    SpanConstants.SPAN_DEFAULT_BACKGROUND_SENSITIVITY
                require(bgSensitivity.isNaN() || 0 < bgSensitivity && bgSensitivity <= 1) {
                    "Illegal background sensitivity: $bgSensitivity, expected range: (0, 1]"
                }
                val clip = if (options.has("clip")) options.valueOf("clip") as Double else SPAN_DEFAULT_CLIP
                require(0 <= clip && clip < 1) { "Illegal clip: $fdr, expected range: [0, 1)" }

                val workingDir = options.valueOf("workdir") as Path
                val id = peaksPath?.stemGz ?: modelPath?.stem ?:
                    reduceIds(listOf(modelId, fdr.toString(), bgSensitivity.toString(), clip.toString()))
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

                if (peaksPath != null) {
                    if (labelsPath != null) {
                        LOG.info("LABELS: $labelsPath")
                        LOG.info("Fdr, background sensitivity, clip options are ignored.")
                    } else {
                        LOG.info("FDR: $fdr")
                        LOG.info("BACKGROUND SENSITIVITY: $bgSensitivity")
                        LOG.info("CLIP: $clip")
                    }
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


                val spanResults = lazySpanResults.value
                val fitInfo = spanResults.fitInfo
                check(fitInfo is AbstractSpanAnalyzeFitInformation) {
                    "Expected ${SpanAnalyzeFitInformation::class.java.simpleName}, got ${fitInfo::class.java.name}"
                }
                val genomeQuery = fitInfo.genomeQuery()
                val fragment = fitInfo.fragment
                val bin = fitInfo.binSize

                if (peaksPath != null) {
                    if (labelsPath == null) {
                        savePeaks(spanResults, fitInfo, genomeQuery, bin, fragment, fdr, bgSensitivity, clip, peaksPath)
                    } else {
                        saveTunedResults(spanResults, genomeQuery, bin, fragment, labelsPath, peaksPath)
                    }
                }
                if (modelPath == null && !keepCacheFiles) {
                    LOG.debug("Clean coverage caches")
                    fitInfo.cleanCaches()
                }
            }
        }
    }


    private fun savePeaks(
        spanResults: SpanFitResults,
        fitInfo: AbstractSpanAnalyzeFitInformation,
        genomeQuery: GenomeQuery,
        bin: Int,
        fragment: Fragment,
        fdr: Double,
        sensitivity: Double,
        clip: Double,
        peaksPath: Path
    ) {
        val peaks = ModelToPeaks.getPeaks(spanResults, genomeQuery, fdr, sensitivity, clip)
        Peak.savePeaks(
            peaks, peaksPath,
            "peak${if (fragment is FixedFragment) "_$fragment" else ""}_${bin}_${fdr}_${sensitivity}_${clip}"
        )
        LOG.info("Saved result to $peaksPath")
        val aboutPeaks = PeaksInfo.compute(
            genomeQuery,
            peaks.map { it.location }.stream(),
            peaksPath.toUri(),
            fitInfo.paths.map { it.treatment }
        )
        val aboutModel = spanResults.modelInformation()
        LOG.info("\n" + (aboutPeaks + aboutModel).joinToString("\n") { (k, v) ->
            "${k.name}: ${k.render(v)}"
        })
    }

    private fun saveTunedResults(
        spanResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        bin: Int,
        fragment: Fragment,
        labelsPath: Path,
        peaksPath: Path
    ) {
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
        val (optimalFDR, optimalSensitivity, optimalClip) = SpanSemiSupervised.PARAMETERS[optimalIndex]
        LOG.info("Optimal settings: FDR=$optimalFDR, SENSITIVITY=$optimalSensitivity, CLIP=$optimalClip")
        labelErrorsGrid.forEachIndexed { i, error ->
            val (fdrTuning, sensitivityTuning, clipTuning) = SpanSemiSupervised.PARAMETERS[i]
            results.addRecord(
                "result",
                "${fdrTuning}_${sensitivityTuning}_${clipTuning}",
                error,
                i == optimalIndex
            )
        }
        results.saveTuningErrors(peaksPath.parent / "${peaksPath.fileName.stem}_errors.csv")
        results.saveOptimalResults(
            peaksPath.parent
                    / "${peaksPath.fileName.stem}_parameters.csv"
        )
        val peaks = ModelToPeaks.getPeaks(
            spanResults, genomeQuery, optimalFDR, optimalSensitivity, optimalClip
        )
        Peak.savePeaks(
            peaks, peaksPath,
            "peak${if (fragment is FixedFragment) "_$fragment" else ""}_" +
                    "${bin}_${optimalFDR}_${optimalSensitivity}"
        )
        LOG.info("Saved result to $peaksPath")
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
    private fun logParametersAndPrepareLazySpanResults(options: OptionSet): Lazy<SpanFitResults> {
        val modelPath = options.valueOf("model") as Path?
        if (modelPath != null) {
            LOG.info("MODEL: $modelPath")
        }
        if (modelPath != null && modelPath.exists && modelPath.size.isNotEmpty()) {
            LOG.debug(
                "Model file $modelPath exists and is not empty, SPAN will use it to substitute " +
                        "the missing command line arguments and verify the provided ones."
            )
            val results = SpanModelFitExperiment.loadResults(modelPath = modelPath)
            check(results.fitInfo is AbstractSpanAnalyzeFitInformation) {
                "Expected ${SpanAnalyzeFitInformation::class.java.simpleName}, got ${results.fitInfo::class.java.name}"
            }
            prepareAndCheckTreatmentControlPaths(options, results.fitInfo, log = true)
            val workingDir = options.valueOf("workdir") as Path
            LOG.info("WORKING DIR: $workingDir")
            val chromSizesPath = options.valueOf("chrom.sizes") as Path?
            if (chromSizesPath != null) {
                checkGenomeInFitInformation(chromSizesPath, results.fitInfo)
            }
            LOG.info("CHROM.SIZES: $chromSizesPath")
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
            // Create fake lazy of already computed results
            return lazyOf(results)
        } else {
            val paths = prepareAndCheckTreatmentControlPaths(options, log = true)
            val workingDir = options.valueOf("workdir") as Path
            LOG.info("WORKING DIR: $workingDir")
            val chromSizesPath = options.valueOf("chrom.sizes") as Path?
            LOG.info("CHROM.SIZES: $chromSizesPath")
            val bin = SpanCLA.getBin(options, log = true)
            val fragment = SpanCLA.getFragment(options, log = true)
            val unique = SpanCLA.getUnique(options, log = true)
            val threshold = SpanCLA.getThreshold(options, log = true)
            val maxIterations = SpanCLA.getMaxIter(options, log = true)
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
            val saveExtendedInfo = options.has("ext")
            LOG.info("EXTENDED MODEL INFO: $saveExtendedInfo")
            val keepCacheFiles = "keep-cache" in options
            LOG.info("KEEP-CACHE: $keepCacheFiles")
            return lazy {
                val genomeQuery = GenomeQuery(Genome[chromSizesPath!!])
                val experiment = getExperimentByModelType(
                    modelType, genomeQuery, paths, unique, fragment, bin,
                    threshold, maxIterations,
                    modelPath, saveExtendedInfo, keepCacheFiles,
                    mapabilityPath
                )
                experiment.results
            }
        }
    }

    private fun getExperimentByModelType(
        modelType: SpanModelType,
        genomeQuery: GenomeQuery,
        paths: List<SpanDataPaths>,
        unique: Boolean,
        fragment: Fragment,
        bin: Int,
        threshold: Double,
        maxIterations: Int,
        modelPath: Path?,
        saveExtendedInfo: Boolean,
        keepCacheFiles: Boolean,
        mapabilityPath: Path?
    ): SpanModelFitExperiment<ClassificationModel, AbstractSpanAnalyzeFitInformation, out Enum<*>> {
        val experiment = when (modelType) {
            SpanModelType.NB2Z_HMM ->
                SpanPeakCallingExperiment.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB2Z_MIXTURE ->
                SpanPeakCallingExperimentNB2ZMixture.getExperiment(
                    genomeQuery, paths, fragment, bin, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB2_HMM ->
                SpanPeakCallingExperimentNB2HMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB3Z_HMM ->
                SpanPeakCallingExperimentNB3ZHMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB5Z_HMM ->
                SpanPeakCallingExperimentNB5ZHMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.NB3_HMM ->
                SpanPeakCallingExperimentNB3HMM.getExperiment(
                    genomeQuery, paths, bin, fragment, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )

            SpanModelType.POISSON_REGRESSION_MIXTURE -> {
                SpanPeakCallingExperimentP2ZRegrMixture.getExperiment(
                    genomeQuery, paths, mapabilityPath, fragment, bin, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
                )
            }

            SpanModelType.NEGBIN_REGRESSION_MIXTURE -> {
                SpanPeakCallingExperimentNB2ZRegrMixture.getExperiment(
                    genomeQuery, paths, mapabilityPath, fragment, bin, unique, modelPath,
                    threshold, maxIterations, saveExtendedInfo, keepCacheFiles
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
