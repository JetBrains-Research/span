package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.experimental.*
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.peaks.PeaksType
import org.jetbrains.bio.span.peaks.getPeaks
import org.jetbrains.bio.span.semisupervised.LocationLabel
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised
import org.jetbrains.bio.span.semisupervised.TuningResults
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.util.*
import org.slf4j.event.Level
import java.nio.file.Path


object SpanCLAAnalyze {

    internal fun analyze(params: Array<String>, experimental: Boolean) {
        with(SpanCLA.getOptionParser()) {
            acceptsAll(
                listOf("t", "treatment"),
                "ChIP-seq treatment file. bam, bed or .bed.gz file;\n" +
                        "If multiple files are given, treated as replicates."
            )
                .requiredUnless("model")
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("c", "control"),
                "Control file. bam, bed or bed.gz file;\n" +
                        "Single control file or separate file per each\n" +
                        "treatment file required."
            )
                .availableIf("treatment")
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.noCheck())
            accepts("labels", "Labels BED file")
                .availableIf("peaks")
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("ms", "multistarts"),
                "Number of multistart runs using different model initialization. Use 0 to disable"
            )
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(Fitter.MULTISTARTS)
            acceptsAll(listOf("msi", "ms-iterations"), "Number of iterations for each multistart run")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(Fitter.MULTISTART_ITERATIONS)
            acceptsAll(listOf("ext"), "Save extended information to model file")
            if (experimental) {
                accepts(
                    "model-type",
                    """
                        Model type.
                        ${SpanModelType.values().joinToString("\n") { "'${it.id}' - ${it.description}" }}
                    """.trimIndent()
                )
                    .withRequiredArg()
                    .defaultsTo(SpanModelType.NB2Z_HMM.id)
                accepts("mapability", "Mapability bigWig file.")
                    .availableIf("treatment")
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.exists())
            }

            parse(params) { options ->
                if ("quiet" in options) {
                    Logs.quiet()
                } else {
                    Logs.addConsoleAppender(if ("debug" in options) Level.DEBUG else Level.INFO)
                }
                SpanCLA.LOG.info("SPAN ${SpanCLA.version()}")
                SpanCLA.LOG.info(
                    "COMMAND: " +
                            "analyze${if (experimental) "-experimental" else ""} ${params.joinToString(" ")}"
                )

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

                val logPath = SpanCLA.configureLogFile(workingDir, id)
                SpanCLA.LOG.info("LOG: $logPath")

                // Call now to preserve params logging order
                val lazySpanResults = peakCallingResults(options, experimental)

                val labelsPath = options.valueOf("labels") as Path?
                val gap = options.valueOf("gap") as Int
                require(gap >= 0) { "Negative gap: $gap" }
                val fdr = options.valueOf("fdr") as Double
                require(0 < fdr && fdr <= 1) { "Illegal fdr: $fdr, expected range: (0, 1)" }
                val peaksTypeCmd = options.valueOf("peaks-type") as String
                require(peaksTypeCmd in PeaksType.values().map { it.cmd }) {
                    "Unknown peaks-type: $peaksTypeCmd, expected ${PeaksType.values().joinToString(", ") { it.cmd }}"
                }
                val peaksType = PeaksType.getTypeFromCmd(peaksTypeCmd)
                val threads = options.valueOf("threads") as Int?
                require(threads == null || threads > 0) {
                    "Not positive threads option: $threads"
                }

                if (peaksPath != null) {
                    if (labelsPath != null) {
                        SpanCLA.LOG.info("LABELS: $labelsPath")
                        SpanCLA.LOG.info("FDR, GAP options are ignored.")
                    } else {
                        SpanCLA.LOG.info("FDR: $fdr")
                        SpanCLA.LOG.info("GAP: $gap")
                    }
                    SpanCLA.LOG.info("TYPE: ${peaksType.cmd}")
                    SpanCLA.LOG.info("PEAKS: $peaksPath")
                } else {
                    SpanCLA.LOG.info("NO peaks path given, process model fitting only.")
                    SpanCLA.LOG.info("LABELS, FDR, GAP options are ignored.")
                }

                configureParallelism(threads)
                SpanCLA.LOG.info("THREADS: ${parallelismLevel()}")

                SpanCLA.checkMemory()

                val spanResults = lazySpanResults.value
                val fitInfo = spanResults.fitInfo
                check(fitInfo is AbstractSpanAnalyzeFitInformation) {
                    "Expected SpanAnalyzeFitInformation, got ${fitInfo::class.java.name}"
                }
                val genomeQuery = fitInfo.genomeQuery()
                val fragment = fitInfo.fragment
                val bin = fitInfo.binSize

                if (peaksPath != null) {
                    if (labelsPath == null) {
                        val peaks = spanResults.getPeaks(genomeQuery, fdr, gap, peaksType = peaksType)
                        Peak.savePeaks(
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
                        SpanCLA.LOG.info("\n" + (aboutPeaks + aboutModel).joinToString("\n") { (k, v) ->
                            "${k.name}: ${k.render(v)}"
                        })
                    } else {
                        val results = TuningResults()
                        SpanCLA.LOG.info("Loading labels $labelsPath...")
                        val labels = LocationLabel.loadLabels(labelsPath, genomeQuery.genome)
                        SpanCLA.LOG.info("Tuning model on the loaded labels...")
                        val (labelErrorsGrid, optimalIndex) = SpanSemiSupervised.tune(
                            spanResults,
                            genomeQuery,
                            labels,
                            "",
                            SpanSemiSupervised.PARAMETERS,
                            peaksType
                        )
                        SpanCLA.LOG.info("Tuning model on the loaded labels complete.")
                        val (optimalFDR, optimalGap) = SpanSemiSupervised.PARAMETERS[optimalIndex]
                        SpanCLA.LOG.info("Optimal settings: FDR=$optimalFDR, GAP=$optimalGap")
                        labelErrorsGrid.forEachIndexed { i, error ->
                            val (fdrTuning, gapTuning) = SpanSemiSupervised.PARAMETERS[i]
                            results.addRecord(
                                "result",
                                "${fdrTuning}_${gapTuning}",
                                error,
                                i == optimalIndex
                            )
                        }
                        results.saveTuningErrors(peaksPath.parent / "${peaksPath.fileName.stem}_errors.csv")
                        results.saveOptimalResults(
                            peaksPath.parent
                                    / "${peaksPath.fileName.stem}_parameters.csv"
                        )
                        val peaks = spanResults.getPeaks(genomeQuery, optimalFDR, optimalGap, peaksType = peaksType)
                        Peak.savePeaks(
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
        options: OptionSet, fitInformation: AbstractSpanAnalyzeFitInformation? = null, log: Boolean = false
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
                "Model file $modelPath exists and is not empty, SPAN will use it to substitute " +
                        "the missing command line arguments and verify the provided ones."
            )
            val results = SpanModelFitExperiment.loadResults(tarPath = modelPath)
            check(results.fitInfo is AbstractSpanAnalyzeFitInformation) {
                "Invalid fit information; expected SpanAnalyzeFitInformation, got ${results.fitInfo::class.java.name}"
            }
            val chromSizesPath = SpanCLA.getAndLogWorkDirAndChromSizes(options, results.fitInfo)
            getAnalyzePaths(options, results.fitInfo, log = true)
            SpanCLA.LOG.info("CHROM.SIZES: $chromSizesPath")
            SpanCLA.getBin(options, results.fitInfo, log = true)
            SpanCLA.getFragment(options, results.fitInfo, log = true)
            SpanCLA.getUnique(options, results.fitInfo, log = true)
            if (results.fitInfo is SpanRegrMixtureAnalyzeFitInformation) {
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
            val threshold = SpanCLA.getThreshold(options, log = true)
            val maxIterations = SpanCLA.getMaxIter(options, log = true)
            val multistarts = getMultistarts(options, log = true)
            val multistartIterations = getMultistartIterations(options, log = true)
            val modelType: SpanModelType
            val mapabilityPath: Path?
            if (experimental) {
                modelType = getModelType(options, modelPath)
                mapabilityPath = if (
                    modelType == SpanModelType.POISSON_REGRESSION_MIXTURE ||
                    modelType == SpanModelType.NEGBIN_REGRESSION_MIXTURE
                ) {
                    getMapabilityPath(options, log = true)
                } else {
                    null
                }
            } else {
                modelType = SpanModelType.NB2Z_HMM
                mapabilityPath = null
            }
            return lazy {
                val saveExtendedInfo = options.has("ext")
                val experiment = when (modelType) {
                    SpanModelType.NB2Z_HMM ->
                        SpanPeakCallingExperiment.getExperiment(
                            genomeQuery, data, bin, fragment, unique, modelPath,
                            threshold, maxIterations, multistarts, multistartIterations,
                            saveExtendedInfo
                        )
                    SpanModelType.NB2Z_MIXTURE ->
                        SpanPeakCallingExperimentNB2ZMixture.getExperiment(
                            genomeQuery, data, fragment, bin, unique, modelPath,
                            threshold, maxIterations, multistarts, multistartIterations,
                            saveExtendedInfo
                        )
                    SpanModelType.NB2_HMM ->
                        SpanPeakCallingExperimentNB2HMM.getExperiment(
                            genomeQuery, data, bin, fragment, unique, modelPath,
                            threshold, maxIterations, multistarts, multistartIterations,
                            saveExtendedInfo
                        )
                    SpanModelType.NB3Z_HMM ->
                        SpanPeakCallingExperimentNB3ZHMM.getExperiment(
                            genomeQuery, data, bin, fragment, unique, modelPath,
                            threshold, maxIterations, multistarts, multistartIterations,
                        )
                    SpanModelType.NB3_HMM3 ->
                        SpanPeakCallingExperimentNB3HMM.getExperiment(
                            genomeQuery, data, bin, fragment, unique, modelPath,
                            threshold, maxIterations, multistarts, multistartIterations,
                        )
                    SpanModelType.POISSON_REGRESSION_MIXTURE -> {
                        SpanPeakCallingExperimentP2ZRegrMixture.getExperiment(
                            genomeQuery, data, mapabilityPath, fragment, bin, unique, modelPath,
                            threshold, maxIterations
                        )
                    }
                    SpanModelType.NEGBIN_REGRESSION_MIXTURE -> {
                        SpanPeakCallingExperimentNB2ZRegrMixture.getExperiment(
                            genomeQuery, data, mapabilityPath, fragment, bin, unique, modelPath,
                            threshold, maxIterations
                        )
                    }
                }
                experiment.results
            }
        }
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
