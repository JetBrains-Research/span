package org.jetbrains.bio.span

import com.google.common.annotations.VisibleForTesting
import joptsimple.OptionParser
import joptsimple.OptionSet
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.FixedFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.experiments.tuning.PeakAnnotation
import org.jetbrains.bio.experiments.tuning.Span
import org.jetbrains.bio.experiments.tuning.TuningResults
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo.compute
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.util.*
import org.jetbrains.bio.util.FileSize.Companion.GB
import java.nio.file.Path

/**
 * Tool for analyzing and comparing ChIP-Seq data.
 * Both procedures rely on the Zero Inflated Negative Binomial Restricted Algorithm.
 *
 * @author Oleg Shpynov
 * @since  14/09/15
 */
@Suppress("UNCHECKED_CAST")
object SpanCLA {
    private val LOG: Logger = Logger.getLogger(SpanCLA::class.java)

    /**
     * Shpynov:
     * Since [Configuration] allows to configure experimentsPath only once,
     * SpanCLA fails to setup correct working directory, if launched within the same process.
     * This is a HACK.
     */
    @VisibleForTesting
    var ignoreConfigurePaths: Boolean = false

    init {
        // Load build properties
        val resource = SpanCLA::class.java.getResource("/span.properties")
        if (resource != null) {
            resource.openStream().use { System.getProperties().load(it) }
        }
    }

    private fun version() =
            "${System.getProperty("span.build.version", "@VERSION@.@build@")} " +
                    "built on ${System.getProperty("span.build.date", "@DATE@")}"


    private const val HELP = """
Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode, experimental
"""
    private const val ANALYZE = "analyze"
    private const val COMPARE = "compare"

    @JvmStatic
    fun main(args: Array<String>) {
        if (args.isEmpty()) {
            System.err.println("ERROR: No command given; $ANALYZE or $COMPARE expected.")
            System.err.println(HELP)
        } else {
            when (args[0]) {
                ANALYZE -> analyze(args.copyOfRange(1, args.size))
                COMPARE -> compare(args.copyOfRange(1, args.size))

                "-?", "-h", "--help" -> println(HELP)
                "-v", "--version" -> println(version())

                else -> {
                    System.err.println("ERROR: Unknown command: ${args[0]}; $ANALYZE or $COMPARE expected.")
                    System.err.println(HELP)
                }
            }
        }
    }


    private fun analyze(params: Array<String>) {
        with(getOptionParser()) {
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
            acceptsAll(listOf("mappability"),
                    "Mappability file. bigWig file;\n" +
                              "treatment file required.")
                    .availableIf("treatment")
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->

                // this value is lazy to ensure the correct logging order
                val lazySpanResults = peakCallingResults(options, params)

                val peaksPath = options.valueOf("peaks") as Path?
                val labelsPath = options.valueOf("labels") as Path?
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val threads = options.valueOf("threads") as Int?

                if (peaksPath != null) {
                    if (labelsPath != null) {
                        LOG.info("LABELS: $labelsPath")
                        LOG.info("FDR, GAP options are ignored.")
                    } else {
                        LOG.info("FDR: $fdr")
                        LOG.info("GAP: $gap")
                    }
                    LOG.info("PEAKS: $peaksPath")
                } else {
                    LOG.info("NO output path given, process model fitting only.")
                    LOG.info("LABELS, FDR, GAP options are ignored.")
                }

                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")

                checkMemory()

                val spanResults = lazySpanResults.value
                val fitInfo = spanResults.fitInfo
                val gq = fitInfo.genomeQuery()
                val fragment = fitInfo.fragment
                val bin = fitInfo.binSize

                if (peaksPath != null) {
                    if (labelsPath == null) {
                        val peaks = spanResults.getPeaks(gq, fdr, gap)
                        savePeaks(
                                peaks, peaksPath,
                                "peak${if (fragment is FixedFragment) "_$fragment" else ""}_${bin}_${fdr}_${gap}"
                        )
                        LOG.info("Saved result to $peaksPath")
                        val aboutPeaks = compute(
                                gq,
                                peaks.map { it.location }.stream(),
                                peaksPath.toUri(),
                                fitInfo.data.map { it.pathTreatment }
                        )
                        val aboutModel = spanResults.about()
                        LOG.info("\n" + (aboutPeaks + aboutModel).map { (k, v) -> "$k: $v" }.joinToString("\n"))
                    } else {
                        val results = TuningResults()
                        LOG.info("Loading labels $labelsPath...")
                        val labels = PeakAnnotation.loadLabels(labelsPath, gq.genome)
                        LOG.info("Tuning model on the loaded labels...")
                        val (labelErrorsGrid, index) = Span.tune(spanResults, labels, "", Span.parameters)
                        LOG.info("Tuning model on the loaded labels complete.")
                        val (optimalFDR, optimalGap) = Span.parameters[index]
                        LOG.info("Optimal settings: FDR=$optimalFDR, GAP=$optimalGap")
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
                        val peaks = spanResults.getPeaks(gq, optimalFDR, optimalGap)
                        savePeaks(
                                peaks, peaksPath,
                                "peak${if (fragment is FixedFragment) "_$fragment" else ""}_" +
                                        "${bin}_${optimalFDR}_$optimalGap"
                        )
                        LOG.info("Saved result to $peaksPath")
                    }
                }
            }
        }
    }

    private fun checkMemory() {
        val maxMemory = Runtime.getRuntime().maxMemory()
        // [Shpynov] This is a hack: since -Xmx4G results in about 3.5G max memory available
        // The reason of such a discrepancy is the size of the garbage collector's survivor space.
        // https://stackoverflow.com/questions/23701207/why-do-xmx-and-runtime-maxmemory-not-agree
        // 3.5 works, however use decreased level to be sure.
        if (maxMemory < 3 * GB) {
            LOG.warn("Recommended memory settings ${FileSize(4 * GB)} are not set. Current settings: ${FileSize(maxMemory)}.\n" +
                    "Please use java memory option '-Xmx4G' to configure memory available for SPAN.")
        }
    }


    private fun compare(params: Array<String>) {
        with(getOptionParser()) {
            acceptsAll(listOf("t1", "treatment1"),
                    "ChIP-seq treatment file 1. bam, bed or .bed.gz file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c1", "control1"),
                    "Control file 1. bam, bed or .bed.gz file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(listOf("t2", "treatment2"),
                    "ChIP-seq treatment file 2. bam, bed or .bed.gz file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c2", "control2"),
                    "Control file 2. bam, bed or .bed.gz file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->

                val workingDir = options.valueOf("workdir") as Path
                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                val treatmentPaths1 = options.valuesOf("treatment1") as List<Path>
                val treatmentPaths2 = options.valuesOf("treatment2") as List<Path>
                val controlPaths1 = options.valuesOf("control1") as List<Path>
                val controlPaths2 = options.valuesOf("control2") as List<Path>
                val mappabilityPaths1 = options.valuesOf("mappability1") as List<Path>
                val mappabilityPaths2 = options.valuesOf("mappability2") as List<Path>
                val fragment = getFragment(options, null, log = false)
                val bin = getBin(options, null, log = false)
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val peaksPath = options.valueOf("peaks") as Path?
                val threads = options.valueOf("threads") as Int?

                // Configure logging
                val id = if (peaksPath != null) {
                    peaksPath.stemGz
                } else {
                    val ids = listOfNotNull(treatmentPaths1, controlPaths1, treatmentPaths2, controlPaths2)
                            .flatMap { paths -> paths.map { it.stemGz } }.toMutableList()
                    ids.add(bin.toString())
                    if (fragment is FixedFragment) {
                        ids.add(fragment.toString())
                    }
                    reduceIds(ids)
                }
                val logPath = configureLogging(
                        "quiet" in options, "debug" in options, id, workingDir
                )


                LOG.info("SPAN ${version()}")
                LOG.info("COMMAND:\ncompare ${params.joinToString(" ")}")
                LOG.info("LOG: $logPath")
                LOG.info("WORKING DIR: $workingDir")
                LOG.info("TREATMENT1: ${treatmentPaths1.joinToString(", ", transform = Path::toString)}")
                if (controlPaths1 != null && controlPaths1.isNotEmpty()) {
                    LOG.info("CONTROL1: ${controlPaths1.joinToString(", ", transform = Path::toString)}")
                } else {
                    LOG.info("CONTROL1: none")
                }
                LOG.info("TREATMENT2: ${treatmentPaths2.joinToString(", ", transform = Path::toString)}")
                if (controlPaths2 != null && controlPaths2.isNotEmpty()) {
                    LOG.info("CONTROL2: ${controlPaths2.joinToString(", ", transform = Path::toString)}")
                } else {
                    LOG.info("CONTROL2: none")
                }
                LOG.info("CHROM.SIZES: $chromSizesPath")
                // Configuration initialization should be configured before any access
                configurePaths(workingDir, chromSizesPath)
                val genome = Genome[chromSizesPath]
                LOG.info("GENOME: ${genome.build}")
                LOG.info("FRAGMENT: $fragment")
                LOG.info("BIN: $bin")
                LOG.info("FDR: $fdr")
                LOG.info("GAP: $gap")
                if (peaksPath != null) {
                    LOG.info("PEAKS: $peaksPath")
                } else {
                    LOG.info("BIN: $bin")
                    LOG.info("NO output path given, process model fitting only.")
                    LOG.info("LABELS, FDR, GAP options are ignored.")
                }
                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")

                checkMemory()


                val gq = GenomeQuery(genome)
                val paths1 = treatmentPaths1.mapIndexed { index, path ->
                    matchTreatmentAndControlsAndMappability(
                            treatmentPaths1[index],
                            controlPaths1[index])
                }
                val coverages1 = treatmentPaths1.map { ReadsQuery(gq, it, fragment = fragment) }
                val paths2 = treatmentPaths2.mapIndexed { index, path ->
                    matchTreatmentAndControlsAndMappability(
                            treatmentPaths2[index],
                            controlPaths2[index])
                }
                val coverages2 = treatmentPaths1.map { ReadsQuery(gq, it, fragment = fragment) }
                val experiment = SpanDifferentialPeakCallingExperiment.getExperiment(
                        gq, paths1, paths2, bin, fragment
                )

                if (peaksPath != null) {
                    val peaks = experiment.results.getPeaks(experiment.genomeQuery, fdr, gap)
                    peaks.forEach { peak ->
                        peak.value =
                                Math.max(1.0, coverages1.map {
                                    it.get().getBothStrandsCoverage(peak.range.on(peak.chromosome))
                                }.average()) /
                                        Math.max(1.0, coverages2.map {
                                            it.get().getBothStrandsCoverage(peak.range.on(peak.chromosome))
                                        }.average())
                    }
                    savePeaks(
                            peaks, peaksPath,
                            "diff${if (fragment is FixedFragment) "_$fragment" else ""}_${bin}_${fdr}_${gap}"
                    )
                    LOG.info("Saved result to $peaksPath")
                } else {
                    experiment.run()
                }
            }
        }
    }


    private fun matchTreatmentAndControlsAndMappability(
            treatmentPaths: Path,
            controlPaths: Path,
            mappabilityPath: Path?
    ): SpanPathsToData {
        return SpanPathsToData(treatmentPaths, controlPaths, mappabilityPath)
    }

    private fun matchTreatmentAndControlsAndMappability(
            treatmentPaths: Path,
            controlPaths: Path
    ): SpanPathsToData {
        return SpanPathsToData(treatmentPaths, controlPaths)
    }

    private fun matchTreatmentAndControls(
            treatmentPaths: List<Path>,
            controlPaths: List<Path>?
    ): List<Pair<Path, Path?>> {
        if (controlPaths != null && controlPaths.isNotEmpty()) {
            if (controlPaths.size != 1 && controlPaths.size != treatmentPaths.size) {
                System.err.println("ERROR: required single control file or separate file for each treatment.")
                System.exit(1)
            }
            return treatmentPaths.zip(Array(treatmentPaths.size) {
                if (controlPaths.size != 1) controlPaths[it] else controlPaths.first()
            })
        }
        return treatmentPaths.map {
            it to null
        }
    }

    private fun configurePaths(outputPath: Path, chromSizesPath: Path?) {
        if (ignoreConfigurePaths) {
            LOG.debug("IGNORE configurePaths")
            return
        }
        outputPath.createDirectories()
        Configuration.experimentsPath = outputPath
        if (chromSizesPath != null) {
            Configuration.genomesPath = chromSizesPath.parent
        }
    }


    private fun getOptionParser(
            bin: Boolean = true,
            fdr: Boolean = true,
            gap: Boolean = true
    ): OptionParser = object : OptionParser() {
        init {
            acceptsAll(listOf("d", "debug"), "Print all the debug information, used for troubleshooting.")
            acceptsAll(listOf("q", "quiet"), "Turn off output")
            acceptsAll(listOf("m", "model"), "Path to model file")
                    .withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                    listOf("cs", "chrom.sizes"),
                    "Chromosome sizes path, can be downloaded at\n" +
                            "http://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes"
            ).requiredUnless("model").withRequiredArg().withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                    listOf("p", "peaks"), "Path to result peaks file in ENCODE broadPeak (BED 6+3) format"
            ).withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                    listOf("fragment"),
                    "Fragment size. If it's an integer, reads are shifted appropriately.\n" +
                            "If it's the string 'auto', the shift is estimated from the data. (default: auto)"
            ).withRequiredArg().withValuesConvertedBy(FragmentConverter())
            if (bin) {
                acceptsAll(listOf("b", "bin"), "Bin size. (default: ${Span.DEFAULT_BIN})")
                        .withRequiredArg()
                        .ofType(Int::class.java)
            }
            if (fdr) {
                acceptsAll(listOf("f", "fdr"), "FDR value.")
                        .availableIf("peaks")
                        .withRequiredArg()
                        .ofType(Double::class.java)
                        .defaultsTo(Span.DEFAULT_FDR)
            }
            if (gap) {
                acceptsAll(listOf("g", "gap"), "Gap size to merge peaks (in bins).")
                        .availableIf("peaks")
                        .withRequiredArg()
                        .ofType(Int::class.java)
                        .defaultsTo(Span.DEFAULT_GAP)
            }
            acceptsAll(listOf("w", "workdir"), "Path to the working dir")
                    .withRequiredArg().withValuesConvertedBy(PathConverter.exists())
                    .defaultsTo(System.getProperty("user.dir").toPath())
            acceptsAll(listOf("threads"), "Parallelism level")
                    .withRequiredArg()
                    .ofType(Int::class.java)
            acceptsAll(listOf("k", "keep-dup"), "Keep duplicates")
                    .withOptionalArg()
                    .ofType(Boolean::class.java)
                    .defaultsTo(true)
        }
    }

    /**
     * Returns path to log file.
     */
    @VisibleForTesting
    internal fun configureLogging(quiet: Boolean, debug: Boolean, id: String, workingDir: Path): Path {
        if (quiet) {
            Logs.quiet()
        }
        // Anyway we configure logging to file.
        Logs.addConsoleAppender(if (debug) Level.DEBUG else Level.INFO)
        val logPath = workingDir / "logs" / "$id.log"
        Logs.addLoggingToFile(logPath)
        return logPath
    }

    /**
     * Configures logging and logs Span version and invocation.
     */
    private fun logGeneralInfoAnalyze(options: OptionSet, params: Array<String>) {
        val peaksPath = options.valueOf("peaks") as Path?
        val modelPath = options.valueOf("model") as Path?
        val workingDir = options.valueOf("workdir") as Path
        // Configure logging
        val id = if (peaksPath != null) {
            peaksPath.stemGz
        } else if (modelPath != null) {
            modelPath.stemGz
        } else {
            // No peaks, no model, generate ID from command-line options.
            // Option parser guarantees that treatment paths are not empty here.
            val (treatmentPaths, controlPaths) = getPaths(options)
            val ids = listOfNotNull(treatmentPaths, controlPaths).flatMap { paths ->
                paths.map { it.stemGz }
            }.toMutableList()
            val bin = getBin(options)
            ids.add(bin.toString())
            val fragment = getFragment(options)
            if (fragment is FixedFragment) {
                ids.add(fragment.size.toString())
            }
            val unique = getUnique(options)
            if (unique) {
                ids.add("unique")
            }
            val labelsPath = options.valueOf("labels") as Path?
            if (labelsPath != null) {
                ids.add(labelsPath.stemGz)
            }
            reduceIds(ids)
        }
        val logPath = configureLogging(
                "quiet" in options, "debug" in options, id, workingDir
        )

        LOG.info("SPAN ${version()}")
        LOG.info("COMMAND:\nanalyze ${params.joinToString(" ")}")
        LOG.info("LOG: $logPath")
    }

    /**
     * Checks supplied command line options against those stored in the fit information.
     * Configures working directory and genomes path (if provided). Logs the progress.
     */
    private fun checkFitInformation(options: OptionSet, fitInformation: SpanFitInformation) {
        getAndLogWorkDirAndChromSizes(options, fitInformation)
        getPaths(options, fitInformation, log = true)
        getBin(options, fitInformation, log = true)
        getFragment(options, fitInformation, log = true)
        getUnique(options, fitInformation, log = true)
    }

    /**
     * Compare the command-line option (null if not given) and stored value (null if not present).
     * Fail if these differ. Return default if both are null. Log the value if [log] is true.
     */
    private fun <T> getProperty(
            commandLineValue: T?, fitInfoValue: T?, default: T, propertyName: String,
            propertyId: String, log: Boolean
    ): T {
        check(fitInfoValue == null || commandLineValue == null || commandLineValue == fitInfoValue) {
            "Stored $propertyName ($fitInfoValue) differs from the command line argument ($commandLineValue)"
        }
        val property = commandLineValue ?: fitInfoValue ?: default
        if (log) {
            LOG.info("$propertyId: $property")
        }
        return property
    }

    private fun getUnique(
            options: OptionSet, fitInformation: SpanFitInformation? = null, log: Boolean = false
    ) = !getProperty(
            if ("keep-dup" in options) options.valueOf("keep-dup") as Boolean else null,
            fitInformation?.unique?.not(), false,
            "'keep duplicates' flag", "KEEP DUPLICATES", log
    )

    private fun getFragment(
            options: OptionSet, fitInformation: SpanFitInformation? = null, log: Boolean = false
    ) = getProperty(
            options.valueOf("fragment") as Fragment?, fitInformation?.fragment, AutoFragment,
            "fragment size", "FRAGMENT", log
    )

    private fun getBin(
            options: OptionSet, fitInformation: SpanFitInformation? = null, log: Boolean = false
    ) = getProperty(
            options.valueOf("bin") as Int?, fitInformation?.binSize, Span.DEFAULT_BIN,
            "bin size", "BIN", log
    )

    private fun getAndLogWorkDirAndChromSizes(
            options: OptionSet, fitInformation: SpanFitInformation? = null
    ): Pair<Path, Path?> {
        val workingDir = options.valueOf("workdir") as Path
        LOG.info("WORKING DIR: $workingDir")

        val chromSizesPath = options.valueOf("chrom.sizes") as Path?
        configurePaths(workingDir, chromSizesPath)

        if (fitInformation != null && chromSizesPath != null) {
            LOG.info("CHROM.SIZES: $chromSizesPath")
            val genome = Genome[chromSizesPath]
            check(genome.build == fitInformation.build) {
                "Stored genome build ${genome.build} differs from " +
                        "the one inferred from the chrom.sizes file $chromSizesPath"
            }
            LOG.info("GENOME: ${genome.build}")
            val chromosomeMap = genome.chromosomeNamesMap
            // we don't check the map equality, since the stored map contains only non-empty chromosomes
            fitInformation.chromosomesSizes.forEach { name, length ->
                check(name in chromosomeMap) {
                    "Stored chromosome $name couldn't be found in the chrom.sizes file $chromSizesPath"
                }
                val chromSizesLength = chromosomeMap.getValue(name).length
                check(chromSizesLength == length) {
                    "Stored chromosome $name length $length differs from $chromSizesLength " +
                            "provided by the chrom.sizes file $chromSizesPath"
                }
            }
        }
        return workingDir to chromSizesPath
    }

    private fun getPaths(
            options: OptionSet, fitInformation: SpanFitInformation? = null, log: Boolean = false
    ): Triple<List<Path>, List<Path>, List<Path>> {
        val commandLineTreatmentPaths = options.valuesOf("treatment") as List<Path>
        val commandLineControlPaths = options.valuesOf("control") as List<Path>
        val commandLineMappabilityPaths = options.valuesOf("mappability") as List<Path>
        if (fitInformation != null && (commandLineTreatmentPaths.isNotEmpty() || commandLineControlPaths.isNotEmpty())) {
            check(commandLineTreatmentPaths.isNotEmpty()) { "Control paths are provided but treatment paths are missing." }
            val paths = commandLineControlPaths.mapIndexed { index, path -> matchTreatmentAndControlsAndMappability(
                    commandLineTreatmentPaths[index],
                    commandLineControlPaths[index],
                    commandLineMappabilityPaths.getOrNull(index))}
            check(paths == fitInformation.data.map { it.pathTreatment to it.pathInput }) {
                "Stored treatment-control pairs ${fitInformation.data.joinToString()} differ from the ones inferred " +
                        "from the command line arguments: ${paths.joinToString()}"
            }

        }
        val treatmentPaths = when {
            fitInformation != null -> fitInformation.data.map { it.pathTreatment }
            commandLineTreatmentPaths.isNotEmpty() -> commandLineTreatmentPaths
            else -> throw IllegalStateException("No treatment files and no existing model file provided, exiting.")
        }
        val controlPaths = when {
            fitInformation != null -> fitInformation.data.mapNotNull { it.pathInput }
            else -> commandLineControlPaths
        }
        if (log) {
            LOG.info("TREATMENT: ${treatmentPaths.joinToString(", ", transform = Path::toString)}")
            if (controlPaths.isNotEmpty()) {
                LOG.info("CONTROL: ${controlPaths.joinToString(", ", transform = Path::toString)}")
            } else {
                LOG.info("CONTROL: none")
            }
        }
        return Triple(treatmentPaths, controlPaths, commandLineMappabilityPaths)
    }

    /**
     * Creates the peak calling experiment using the supplied command line arguments.
     * Logs the progress.
     */
    private fun constructPeakCallingExperiment(
            options: OptionSet
    ): Lazy<SpanPeakCallingExperiment<out ClassificationModel, ZLH>> {
        val (_, chromSizesPath) = getAndLogWorkDirAndChromSizes(options)
        // option parser guarantees that chrom.sizes are not null here
        val genomeQuery = GenomeQuery(Genome[chromSizesPath!!])
        val (treatmentPaths, controlPaths, mappabilityPaths) = getPaths(options, log = true)
        val data = controlPaths.mapIndexed { index, path ->
            matchTreatmentAndControlsAndMappability(
                    treatmentPaths[index],
                    controlPaths[index],
                    mappabilityPaths.getOrNull(index))
        }
        val bin = getBin(options, log = true)
        val fragment = getFragment(options, log = true)
        val unique = getUnique(options, log = true)
        val modelPath = options.valueOf("model") as Path?
        return lazy { SpanPeakCallingExperiment.getExperiment(genomeQuery, data, bin, fragment, unique, modelPath) }
    }

    /**
     * Configure logging and get [SpanFitResults] in a most concise and effective way.
     * Parses and logs most of the command line arguments.
     * Doesn't fit the model; this happens only if .value is invoked
     */
    private fun peakCallingResults(options: OptionSet, params: Array<String>): Lazy<SpanFitResults> {
        logGeneralInfoAnalyze(options, params)

        val modelPath = options.valueOf("model") as Path?
        if (modelPath != null) {
            LOG.info("MODEL: $modelPath")
            if (modelPath.exists && modelPath.size.isNotEmpty()) {
                LOG.debug(
                        "Model file $modelPath exists and is not empty, Span will use it to substitute " +
                                "the missing command line arguments and verify the provided ones."
                )
                val results = SpanModelFitExperiment.loadResults(tarPath = modelPath)
                checkFitInformation(options, results.fitInfo)
                return lazyOf(results)
            }
        }
        val experiment = constructPeakCallingExperiment(options)
        return lazy { experiment.value.results }
    }
}