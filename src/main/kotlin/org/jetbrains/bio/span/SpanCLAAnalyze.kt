package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.experiments.tuning.LocationLabel
import org.jetbrains.bio.experiments.tuning.TuningResults
import org.jetbrains.bio.experiments.tuning.tools.Span
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.genome.query.stemGz
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.util.*
import java.nio.file.Path


object SpanCLAAnalyze {

    internal fun analyze(params: Array<String>, experimental: Boolean) {
        with(SpanCLA.getOptionParser()) {
            acceptsAll(listOf("t", "treatment"),
                    "ChIP-seq treatment file. bam, bed or .bed.gz file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .requiredUnless("model")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(listOf("c", "control"),
                    "Control file. bam, bed or bed.gz file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .availableIf("treatment")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(listOf("labels"), "Labels BED file")
                    .availableIf("peaks")
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("ms", "multistarts"),
                    "Number of multistart runs using different model initialization. Use 0 to disable")
                    .withRequiredArg()
                    .ofType(Int::class.java)
                    .defaultsTo(Fitter.MULTISTARTS)
            acceptsAll(listOf("msi", "ms-iterations"), "Number of iterations for each multistart run")
                    .withRequiredArg()
                    .ofType(Int::class.java)
                    .defaultsTo(Fitter.MULTISTART_ITERATIONS)

            if (experimental) {
                acceptsAll(
                        listOf("type"),
                        "Model type. Either 'nbhmm' for negative binomial HMM (default), " +
                                "'prm' for Poisson regression mixture (experimental) " +
                                "or 'nbrm' for Negative Binomial regression mixture (experimental)."
                ).withRequiredArg()
                acceptsAll(listOf("mapability"), "Mapability bigWig file.")
                        .availableIf("treatment")
                        .withRequiredArg()
                        .withValuesConvertedBy(PathConverter.exists())
            }

            SpanCLA.LOG.info("SPAN ${SpanCLA.version()}")
            SpanCLA.LOG.info("COMMAND:\nanalyze ${params.joinToString(" ")}")
            parse(params) { options ->

                val peaksPath = options.valueOf("peaks") as Path?
                val modelPath = options.valueOf("model") as Path?
                val workingDir = options.valueOf("workdir") as Path
                // Configure logging
                val id = peaksPath?.stemGz ?: if (modelPath != null) {
                    modelPath.stem
                } else {
                    // No peaks, no model, generate ID from command-line options.
                    // Option parser guarantees that treatment paths are not empty here.
                    val data = getAnalyzePaths(options)
                    val ids = listOf(data.map { it.treatment }, data.mapNotNull { it.control }).flatMap { paths ->
                        paths.map { it.stemGz }
                    }.toMutableList()
                    val bin = SpanCLA.getBin(options)
                    ids.add(bin.toString())
                    val fragment = SpanCLA.getFragment(options)
                    if (fragment is FixedFragment) {
                        ids.add(fragment.size.toString())
                    }
                    val unique = SpanCLA.getUnique(options)
                    if (unique) {
                        ids.add("unique")
                    }
                    val labelsPath = options.valueOf("labels") as Path?
                    if (labelsPath != null) {
                        ids.add(labelsPath.stemGz)
                    }
                    reduceIds(ids)
                }

                val logPath = SpanCLA.configureLogging("quiet" in options, "debug" in options, id, workingDir)
                SpanCLA.LOG.info("LOG: $logPath")

                // Call now to preserve params logging order
                val lazySpanResults = peakCallingResults(options, experimental)

                val labelsPath = options.valueOf("labels") as Path?
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val threads = options.valueOf("threads") as Int?

                if (peaksPath != null) {
                    if (labelsPath != null) {
                        SpanCLA.LOG.info("LABELS: $labelsPath")
                        SpanCLA.LOG.info("FDR, GAP options are ignored.")
                    } else {
                        SpanCLA.LOG.info("FDR: $fdr")
                        SpanCLA.LOG.info("GAP: $gap")
                    }
                    SpanCLA.LOG.info("PEAKS: $peaksPath")
                } else {
                    SpanCLA.LOG.info("NO output path given, process model fitting only.")
                    SpanCLA.LOG.info("LABELS, FDR, GAP options are ignored.")
                }

                configureParallelism(threads)
                SpanCLA.LOG.info("THREADS: ${parallelismLevel()}")

                SpanCLA.checkMemory()

                val spanResults = lazySpanResults.value
                val fitInfo = spanResults.fitInfo
                check(fitInfo is SpanAnalyzeFitInformation) {
                    "Expected SpanAnalyzeFitInformation, got ${fitInfo::class.java.name}"
                }
                val genomeQuery = fitInfo.genomeQuery()
                val fragment = fitInfo.fragment
                val bin = fitInfo.binSize

                if (peaksPath != null) {
                    if (labelsPath == null) {
                        val peaks = spanResults.getPeaks(genomeQuery, fdr, gap)
                        savePeaks(
                                peaks, peaksPath,
                                "peak${if (fragment is FixedFragment) "_$fragment" else ""}_${bin}_${fdr}_${gap}"
                        )
                        SpanCLA.LOG.info("Saved result to $peaksPath")
                        val aboutPeaks = PeaksInfo.compute(
                                genomeQuery,
                                peaks.map { it.location }.stream(),
                                peaksPath.toUri(),
                                fitInfo.data.map { it.treatment }
                        )
                        val aboutModel = spanResults.about()
                        SpanCLA.LOG.info("\n" + (aboutPeaks + aboutModel).map { (k, v) ->
                            "${k.name}: ${k.render(v)}"
                        }.joinToString("\n"))
                    } else {
                        val results = TuningResults()
                        SpanCLA.LOG.info("Loading labels $labelsPath...")
                        val labels = LocationLabel.loadLabels(labelsPath, genomeQuery.genome)
                        SpanCLA.LOG.info("Tuning model on the loaded labels...")
                        val (labelErrorsGrid, index) = Span.tune(spanResults, labels, "", Span.parameters)
                        SpanCLA.LOG.info("Tuning model on the loaded labels complete.")
                        val (optimalFDR, optimalGap) = Span.parameters[index]
                        SpanCLA.LOG.info("Optimal settings: FDR=$optimalFDR, GAP=$optimalGap")
                        labelErrorsGrid.forEachIndexed { i, error ->
                            results.addRecord(
                                    "result",
                                    Span.transform(Span.parameters[i]),
                                    error,
                                    i == index
                            )
                        }
                        results.saveTuningErrors(peaksPath.parent / "${peaksPath.fileName.stem}_errors.csv")
                        results.saveOptimalResults(peaksPath.parent
                                / "${peaksPath.fileName.stem}_parameters.csv")
                        val peaks = spanResults.getPeaks(genomeQuery, optimalFDR, optimalGap)
                        savePeaks(
                                peaks, peaksPath,
                                "peak${if (fragment is FixedFragment) "_$fragment" else ""}_" +
                                        "${bin}_${optimalFDR}_$optimalGap"
                        )
                        SpanCLA.LOG.info("Saved result to $peaksPath")
                    }
                }
            }
        }
    }

    /**
     * Retrieves the paths (treatment, optional control, and optional mapability)
     * either from command-line options or from the stored fit information.
     *
     * If both are available, checks that they are consistent.
     */
    private fun getAnalyzePaths(
            options: OptionSet, fitInformation: SpanAnalyzeFitInformation? = null, log: Boolean = false
    ): List<SpanDataPaths> {
        val commandLineTreatmentPaths = options.valuesOf("treatment") as List<Path>
        val commandLineControlPaths = options.valuesOf("control") as List<Path>

        val commandLinePaths = SpanCLA.getCommandLinePaths(
                commandLineTreatmentPaths, commandLineControlPaths
        )

        val fitInfoPaths = fitInformation?.data
        if (commandLinePaths != null && fitInfoPaths != null && commandLinePaths != fitInfoPaths) {
            throw IllegalStateException(
                    "Stored treatment-control pairs ${fitInfoPaths.joinToString()} differ from the ones inferred " +
                            "from the command line arguments: ${commandLinePaths.joinToString()}"
            )
        }

        val paths = commandLinePaths
                ?: fitInfoPaths
                ?: throw IllegalStateException("No treatment files and no existing model file provided, exiting.")

        if (log) {
            SpanCLA.LOG.info("TREATMENT: ${paths.map { it.treatment }.joinToString(", ", transform = Path::toString)}")
            paths.mapNotNull { it.control }.let {
                if (it.isNotEmpty()) {
                    SpanCLA.LOG.info("CONTROL: ${it.joinToString(", ", transform = Path::toString)}")
                } else {
                    SpanCLA.LOG.info("CONTROL: none")
                }
            }
        }
        return paths
    }

    internal fun getMultistarts(
            options: OptionSet, log: Boolean = false
    ) = SpanCLA.getProperty(
            options.valueOf("multistarts") as Int?, null, Fitter.MULTISTARTS,
            "multistarts", "MULTISTARTS", log
    )

    internal fun getMultistartIterations(
            options: OptionSet, log: Boolean = false
    ) = SpanCLA.getProperty(
            options.valueOf("ms-iterations") as Int?, null, Fitter.MULTISTART_ITERATIONS,
            "multistart iterations", "MULTISTART ITERATIONS", log
    )

    /**
     * Configure logging and get [SpanFitResults] in a most concise and effective way.
     * Parses and logs most of the command line arguments.
     * Doesn't fit the model; this happens only if .value is invoked
     */
    private fun peakCallingResults(options: OptionSet, experimental: Boolean): Lazy<SpanFitResults> {
        val modelPath = options.valueOf("model") as Path?
        if (modelPath != null) {
            require(modelPath.extension == "span") {
                "Unrecognized model extension '.${modelPath.extension}', should be '.span'."
            }
            SpanCLA.LOG.info("MODEL: $modelPath")
        }
        if (modelPath != null && modelPath.exists && modelPath.size.isNotEmpty()) {
            SpanCLA.LOG.debug(
                    "Model file $modelPath exists and is not empty, Span will use it to substitute " +
                            "the missing command line arguments and verify the provided ones."
            )
            val results = SpanModelFitExperiment.loadResults(tarPath = modelPath)
            check(results.fitInfo is SpanAnalyzeFitInformation) {
                "Invalid fit information; expected SpanAnalyzeFitInformation, got ${results.fitInfo::class.java.name}"
            }
            val chromSizesPath = SpanCLA.getAndLogWorkDirAndChromSizes(options, results.fitInfo)
            getAnalyzePaths(options, results.fitInfo, log = true)
            SpanCLA.LOG.info("CHROM.SIZES: $chromSizesPath")
            SpanCLA.getBin(options, results.fitInfo, log = true)
            SpanCLA.getFragment(options, results.fitInfo, log = true)
            SpanCLA.getUnique(options, results.fitInfo, log = true)
            if (results.fitInfo is Span2FitInformation) {
                getMapabilityPath(options, results.fitInfo, log = true)
            }
            return lazyOf(results)
        } else {
            val chromSizesPath = SpanCLA.getAndLogWorkDirAndChromSizes(options)
            val genomeQuery = GenomeQuery(Genome[chromSizesPath!!])
            val data = getAnalyzePaths(options, log = true)
            SpanCLA.LOG.info("CHROM.SIZES: $chromSizesPath")
            val bin = SpanCLA.getBin(options, log = true)
            val fragment = SpanCLA.getFragment(options, log = true)
            val unique = SpanCLA.getUnique(options, log = true)
            val maxIter = SpanCLA.getMaxIter(options, log = true)
            val multistarts = getMultistarts(options, log = true)
            val multistartIterations = getMultistartIterations(options, log = true)
            val threshold = SpanCLA.getThreshold(options, log = true)
            val modelType: SpanModel
            val mapabilityPath: Path?
            if (experimental) {
                modelType = getModelType(options, modelPath)
                mapabilityPath = if (
                        modelType == SpanModel.POISSON_REGRESSION_MIXTURE ||
                        modelType == SpanModel.NEGBIN_REGRESSION_MIXTURE
                ) {
                    getMapabilityPath(options, log = true)
                } else {
                    null
                }
            } else {
                modelType = SpanModel.NB_HMM
                mapabilityPath = null
            }
            return lazy {
                val experiment = when (modelType) {
                    SpanModel.NB_HMM ->
                        SpanPeakCallingExperiment.getExperiment(
                                genomeQuery, data, bin, fragment, unique, modelPath,
                                threshold, maxIter, multistarts, multistartIterations
                        )
                    SpanModel.POISSON_REGRESSION_MIXTURE -> {
                        if (multistarts > 0) {
                            SpanCLA.LOG.error("Multistart is not supported for $modelType")
                        }
                        Span2PeakCallingExperiment.getExperiment(
                                genomeQuery, data, mapabilityPath, fragment, bin, unique, modelPath,
                                threshold, maxIter
                        )
                    }
                    SpanModel.NEGBIN_REGRESSION_MIXTURE -> {
                        if (multistarts > 0) {
                            SpanCLA.LOG.error("Multistart is not supported for $modelType")
                        }
                        Span3PeakCallingExperiment.getExperiment(
                                genomeQuery, data, mapabilityPath, fragment, bin, unique, modelPath,
                                threshold, maxIter
                        )
                    }

                }
                experiment.results
            }
        }
    }

    private fun getMapabilityPath(
            options: OptionSet, fitInfo: Span2FitInformation? = null, log: Boolean = false
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
            options.valueOf("type")?.let {
                when (it) {
                    "nbhmm" -> SpanModel.NB_HMM
                    "prm" -> SpanModel.POISSON_REGRESSION_MIXTURE
                    "nbrm" -> SpanModel.NEGBIN_REGRESSION_MIXTURE
                    else -> throw IllegalArgumentException("Unrecognized value for --type command line option: $it")
                }
            },
            modelPath?.let {
                when (it.extension) {
                    "span" -> SpanModel.NB_HMM
                    "span2" -> SpanModel.POISSON_REGRESSION_MIXTURE
                    "span3" -> SpanModel.NEGBIN_REGRESSION_MIXTURE
                    else -> throw IllegalArgumentException(
                            "Unrecognized model extension '.${it.extension}', should be either '.span' or '.span2'."
                    )
                }
            },
            SpanModel.NB_HMM,
            "model type", "MODEL TYPE", true
    )

}

enum class SpanModel(val description: String) {
    NB_HMM("negative binomial HMM"),
    POISSON_REGRESSION_MIXTURE("Poisson regression mixture"),
    NEGBIN_REGRESSION_MIXTURE("Negative Binomial Regression mixture");

    override fun toString() = description
}
