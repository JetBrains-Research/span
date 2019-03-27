package org.jetbrains.bio.span

import com.google.common.annotations.VisibleForTesting
import joptsimple.OptionParser
import org.apache.log4j.Level
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.coverage.AutoFragment
import org.jetbrains.bio.coverage.FixedFragment
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.experiments.fit.SpanDifferentialPeakCallingExperiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.experiments.tuning.PeakAnnotation
import org.jetbrains.bio.experiments.tuning.Span
import org.jetbrains.bio.experiments.tuning.TuningResults
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.PeaksInfo.compute
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
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
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c", "control"),
                    "Control file. bam, bed or bed.gz file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("labels"), "Labels BED file")
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->

                val workingDir = options.valueOf("workdir") as Path
                val treatmentPaths = (options.valuesOf("treatment") as List<Path>)
                val controlPaths = (options.valuesOf("control") as List<Path>?)
                val peaksPath = options.valueOf("peaks") as Path?
                val labelsPath = options.valueOf("labels") as Path?
                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                val bin = options.valueOf("bin") as Int
                val fragment = options.valueOf("fragment") as Fragment
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val threads = options.valueOf("threads") as Int?
                val unique = "keep-dup" !in options
                val modelPath = if ("model" in options) options.valueOf("model") as Path else null

                // Configure logging
                val id = if (peaksPath != null) {
                    peaksPath.stemGz
                } else {
                    val ids = listOfNotNull(treatmentPaths, controlPaths).flatMap { paths ->
                        paths.map { it.stemGz }
                    }.toMutableList()
                    ids.add(bin.toString())
                    if (fragment is FixedFragment) {
                        ids.add(fragment.size.toString())
                    }
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
                LOG.info("WORKING DIR: $workingDir")
                LOG.info("TREATMENT: ${treatmentPaths.joinToString(", ", transform = Path::toString)}")
                if (controlPaths != null && controlPaths.isNotEmpty()) {
                    LOG.info("CONTROL: ${controlPaths.joinToString(", ", transform = Path::toString)}")
                } else {
                    LOG.info("CONTROL: none")
                }
                // Configuration initialization should be configured before any access
                configurePaths(workingDir, chromSizesPath)
                val genome = Genome[chromSizesPath]
                LOG.info("CHROM.SIZES: $chromSizesPath")
                LOG.info("GENOME: ${genome.build}")
                if (modelPath != null) {
                    LOG.info("MODEL: $modelPath")
                }
                LOG.info("FRAGMENT: $fragment")
                if (!unique) {
                    LOG.info("KEEP DUPLICATES: true")
                }
                LOG.info("BIN: $bin")
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

                val gq = GenomeQuery(genome)
                val paths = matchTreatmentAndControls(treatmentPaths, controlPaths)
                val peakCallingExperiment = SpanPeakCallingExperiment.getExperiment(
                    gq, paths, bin, fragment, unique, modelPath
                )
                if (peaksPath != null) {
                    if (labelsPath == null) {
                        val peaks = peakCallingExperiment.results.getPeaks(gq, fdr, gap)
                        savePeaks(
                            peaks, peaksPath,
                            "peak${if (fragment is FixedFragment) "_$fragment" else ""}_${bin}_${fdr}_${gap}"
                        )
                        LOG.info("Saved result to $peaksPath")
                        val aboutPeaks = compute(
                            gq,
                                peaks.map { it.location }.stream(),
                                peaksPath.toUri(),
                                peakCallingExperiment.fitInformation.data.map { it.path.toPath() }
                        )
                        val aboutModel = peakCallingExperiment.results.about()
                        LOG.info("\n" + (aboutPeaks + aboutModel).map { (k, v) -> "$k: $v" }.joinToString("\n"))
                    } else {
                        val results = TuningResults()
                        val labels = PeakAnnotation.loadLabels(labelsPath, gq.genome)
                        val (labelErrorsGrid, index) = Span.tune(peakCallingExperiment, labels, "", Span.parameters)
                        val (optimalFDR, optimalGap) = Span.parameters[index]
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
                        val peaks = peakCallingExperiment.results.getPeaks(gq, optimalFDR, optimalGap)
                        savePeaks(
                            peaks, peaksPath,
                            "peak${if (fragment is FixedFragment) "_$fragment" else ""}_" +
                                "${bin}_${optimalFDR}_$optimalGap"
                        )
                        LOG.info("Saved result to $peaksPath")
                    }
                } else {
                    // Just fit the model
                    peakCallingExperiment.run()
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
                val controlPaths1 = options.valuesOf("control1") as List<Path>?
                val controlPaths2 = options.valuesOf("control2") as List<Path>?
                val fragment = options.valueOf("fragment") as Fragment
                val bin = options.valueOf("bin") as Int
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
                val paths1 = matchTreatmentAndControls(treatmentPaths1, controlPaths1)
                val coverages1 = treatmentPaths1.map { ReadsQuery(gq, it, fragment = fragment) }
                val paths2 = matchTreatmentAndControls(treatmentPaths2, controlPaths2)
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

    private fun configurePaths(outputPath: Path, chromSizesPath: Path) {
        if (ignoreConfigurePaths) {
            LOG.debug("IGNORE configurePaths")
            return
        }
        outputPath.createDirectories()
        Configuration.experimentsPath = outputPath
        Configuration.genomesPath = chromSizesPath.parent
    }


    private fun getOptionParser(
            bin: Boolean = true,
            fdr: Boolean = true,
            gap: Boolean = true
    ): OptionParser = object : OptionParser() {
        init {
            acceptsAll(listOf("d", "debug"), "Print all the debug information, used for troubleshooting.")
            acceptsAll(listOf("q", "quiet"), "Turn off output")
            acceptsAll(
                listOf("cs", "chrom.sizes"),
                "Chromosome sizes path, can be downloaded at\n" +
                    "http://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes"
            ).withRequiredArg().required().withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("p", "peaks"), "Path to result peaks file in ENCODE broadPeak (BED 6+3) format"
            ).withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("fragment"),
                "Fragment size. If it's an integer, reads are shifted appropriately.\n" +
                        "If it's the string 'auto', the shift is estimated from the data."
            ).withRequiredArg().withValuesConvertedBy(FragmentConverter()).defaultsTo(AutoFragment)
            if (bin) {
                acceptsAll(listOf("b", "bin"), "Bin size.")
                        .withRequiredArg()
                        .ofType(Int::class.java)
                        .defaultsTo(Span.DEFAULT_BIN)
            }
            if (fdr) {
                acceptsAll(listOf("f", "fdr"), "FDR value.")
                        .withRequiredArg()
                        .ofType(Double::class.java)
                        .defaultsTo(Span.DEFAULT_FDR)
            }
            if (gap) {
                acceptsAll(listOf("g", "gap"), "Gap size to merge peaks (in bins).")
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
            acceptsAll(listOf("m", "model"), "Path to model file")
                    .withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
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
}
