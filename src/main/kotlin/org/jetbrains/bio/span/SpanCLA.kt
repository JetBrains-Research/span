package org.jetbrains.bio.span

import com.google.common.annotations.VisibleForTesting
import com.google.common.io.ByteStreams
import joptsimple.BuiltinHelpFormatter
import joptsimple.OptionParser
import org.apache.log4j.Level
import org.apache.log4j.LogManager
import org.apache.log4j.Logger
import org.jetbrains.bio.Configuration
import org.jetbrains.bio.Logs
import org.jetbrains.bio.experiments.histones.PeaksInfo
import org.jetbrains.bio.experiments.tuning.PeakAnnotation
import org.jetbrains.bio.experiments.tuning.SPAN
import org.jetbrains.bio.experiments.tuning.TuningResults
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.BinnedReadsQuery
import org.jetbrains.bio.query.readsName
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.util.*
import org.jetbrains.bio.util.FileSize.Companion.GB
import java.io.PrintStream
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
    private val LOG: Logger

    /**
     * Shpynov:
     * Since [Configuration] allows to configure experimentsPath only once,
     * SpanCLA fails to setup correct working directory, if launched within the same process.
     * This is a HACK.
     */
    @VisibleForTesting
    internal var ignoreConfigurePaths: Boolean = false

    init {
        // Add appender before initializing logger to avoid Log4j warnings
        Logs.addConsoleAppender(Level.INFO)
        LOG = Logger.getLogger(SpanCLA::class.java)

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
                    "ChIP-seq treatment file. bam, bed, .bed.gz or bigWig file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c", "control"),
                    "Control file. bam, bed, bed.gz or bigWig file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("labels"), "Labels BED file")
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->
                configureLogging("quiet" in options, "debug" in options)

                val workingDir = options.valueOf("workdir") as Path
                // We would like to reuse as much caching as possible, resolve all the symbolic links
                val treatmentPaths = (options.valuesOf("treatment") as List<Path>).map { it.toRealPath() }
                val controlPaths = (options.valuesOf("control") as List<Path>?)?.map { it.toRealPath() }
                val outputBed = options.valueOf("output") as Path?
                val labelsPath = options.valueOf("labels") as Path?
                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                val chromosomes = options.valuesOf("only") as List<String>
                val bin = options.valueOf("bin") as Int
                val fragment = options.valueOf("fragment") as Int?
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val threads = options.valueOf("threads") as Int?

                // Configure logs path
                val logPath: Path
                if (outputBed != null) {
                    logPath = workingDir / "logs" / "${outputBed.readsName()}.log"
                } else {
                    var ids = listOfNotNull(treatmentPaths, controlPaths).flatMap { it.map { it.readsName() } }
                    ids += bin.toString()
                    if (labelsPath != null) {
                        ids += labelsPath.readsName()
                    }
                    logPath = workingDir / "logs"/ "${reduceIds(ids)}.log"
                }
                Logs.addLoggingToFile(logPath)

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
                val genomeQuery = loadGenomeQuery(chromSizesPath, chromosomes)
                LOG.info("CHROM.SIZES: $chromSizesPath")
                LOG.info("GENOME: ${genomeQuery.id}")
                LOG.info("FRAGMENT: $fragment")
                if (outputBed != null) {
                    if (labelsPath != null) {
                        LOG.info("LABELS: $labelsPath")
                        LOG.info("BIN, FDR, GAP options are ignored.")
                    } else {
                        LOG.info("BIN: $bin")
                        LOG.info("FDR: $fdr")
                        LOG.info("GAP: $gap")
                    }
                    LOG.info("OUTPUT: $outputBed")
                } else {
                    LOG.info("BIN: $bin")
                    LOG.info("NO output path given, process model fitting only.")
                    LOG.info("LABELS, FDR, GAP options are ignored.")
                }

                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")

                checkMemory()

                configurePaths(workingDir, chromSizesPath)
                val coverageQueries = coverageQueries(genomeQuery, treatmentPaths, controlPaths, bin, fragment)

                val peakCallingExperiment = Span.getPeakCallingExperiment(genomeQuery, coverageQueries, bin)
                if (outputBed != null) {
                    if (labelsPath == null) {
                        val peaks = peakCallingExperiment.results.getPeaks(genomeQuery, fdr, gap)
                        savePeaks(peaks, outputBed,
                                "peak${if (fragment != null) "_$fragment" else ""}_${bin}_${fdr}_${gap}")
                        LOG.info("Saved result to $outputBed")
                        LOG.info("\n" + PeaksInfo.aboutText(genomeQuery,
                                peaks.map { it.location }.stream(),
                                outputBed.toUri(),
                                coverageQueries))
                    } else {
                        val results = TuningResults()
                        val labels = PeakAnnotation.loadLabels(labelsPath, genomeQuery.build)
                        val (labelErrorsGrid, index) = SPAN.tune(peakCallingExperiment, labels, "", SPAN.parameters)
                        val (optimalFDR, optimalGap) = SPAN.parameters[index]
                        labelErrorsGrid.forEachIndexed { i, error ->
                            results.addRecord("result",
                                    SPAN.transform(SPAN.parameters[i]),
                                    error,
                                    i == index)
                        }
                        results.saveTuningErrors(outputBed.parent / "${outputBed.fileName.stem}_errors.csv")
                        results.saveOptimalResults(outputBed.parent
                                / "${outputBed.fileName.stem}_parameters.csv")
                        val peaks = peakCallingExperiment.results.getPeaks(genomeQuery, optimalFDR, optimalGap)
                        savePeaks(peaks, outputBed, "peak${if (fragment != null) "_$fragment" else ""}_" +
                                "${bin}_${optimalFDR}_$optimalGap")
                        LOG.info("Saved result to $outputBed")
                        LOG.info("\n" + PeaksInfo.aboutText(genomeQuery,
                                peaks.map { it.location }.stream(),
                                outputBed.toUri(),
                                coverageQueries))
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
                    "ChIP-seq treatment file 1. bam, bed, .bed.gz or bigWig file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c1", "control1"),
                    "Control file 1. bam, bed, .bed.gz or bigWig file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(listOf("t2", "treatment2"),
                    "ChIP-seq treatment file 2. bam, bed, .bed.gz or bigWig file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c2", "control2"),
                    "Control file 2. bam, bed, .bed.gz or bigWig file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->
                configureLogging("quiet" in options, "debug" in options)

                val workingDir = options.valueOf("workdir") as Path
                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                val chromosomes = options.valuesOf("only") as List<String>
                val treatmentPaths1 = options.valuesOf("treatment1") as List<Path>
                val treatmentPaths2 = options.valuesOf("treatment2") as List<Path>
                val controlPaths1 = options.valuesOf("control1") as List<Path>?
                val controlPaths2 = options.valuesOf("control2") as List<Path>?
                val fragment = options.valueOf("fragment") as Int?
                val bin = options.valueOf("bin") as Int
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val outputBed = options.valueOf("output") as Path?
                val threads = options.valueOf("threads") as Int?

                // Configure logs path
                val logPath: Path
                if (outputBed != null) {
                    logPath = workingDir / "logs" / "${outputBed.readsName()}.log"
                } else {
                    var ids = listOfNotNull(treatmentPaths1, controlPaths1, treatmentPaths2, controlPaths2)
                            .flatMap { it.map { it.readsName() } }
                    ids += bin.toString()
                    logPath = workingDir / "logs" / "${reduceIds(ids)}.log"
                }
                Logs.addLoggingToFile(logPath)

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
                val genomeQuery = loadGenomeQuery(chromSizesPath, chromosomes)
                LOG.info("GENOME: ${genomeQuery.id}")
                LOG.info("FRAGMENT: $fragment")
                LOG.info("BIN: $bin")
                LOG.info("FDR: $fdr")
                LOG.info("GAP: $gap")
                if (outputBed != null) {
                    LOG.info("OUTPUT: $outputBed")
                } else {
                    LOG.info("BIN: $bin")
                    LOG.info("NO output path given, process model fitting only.")
                    LOG.info("LABELS, FDR, GAP options are ignored.")
                }
                configureParallelism(threads)
                LOG.info("THREADS: ${parallelismLevel()}")

                checkMemory()

                configurePaths(workingDir, chromSizesPath)

                val coverageQueries1 = coverageQueries(genomeQuery, treatmentPaths1, controlPaths1, bin, fragment)
                val coverageQueries2 = coverageQueries(genomeQuery, treatmentPaths2, controlPaths2, bin, fragment)
                val experiment = Span.getDifferentialPeakCallingExperiment(
                        genomeQuery, coverageQueries1, coverageQueries2, bin)
                if (outputBed != null) {
                    val peaks = experiment.results.getPeaks(experiment.genomeQuery, fdr, gap)
                    peaks.forEach { peak ->
                        peak.value =
                                Math.max(1.0, coverageQueries1.map {
                                    it.get().getBothStrandCoverage(peak.range.on(peak.chromosome))
                                }.average()) /
                                Math.max(1.0, coverageQueries2.map {
                                    it.get().getBothStrandCoverage(peak.range.on(peak.chromosome))
                                }.average())
                    }
                    savePeaks(peaks, outputBed, "diff${if (fragment != null) "_$fragment" else ""}_${bin}_${fdr}_${gap}")
                    LOG.info("Saved result to $outputBed")
                } else {
                    experiment.run()
                }
            }
        }
    }


    private fun coverageQueries(genomeQuery: GenomeQuery,
                                treatmentPaths: List<Path>,
                                controlPaths: List<Path>?,
                                bin: Int,
                                fragment: Int?): List<BinnedReadsQuery> {
        if (controlPaths != null && controlPaths.isNotEmpty()) {
            return if (controlPaths.size == 1) {
                treatmentPaths.map {
                    BinnedReadsQuery(genomeQuery, it, bin, controlPaths.first(), unique = true, fragment = fragment)
                }
            } else {
                if (controlPaths.size != treatmentPaths.size) {
                    System.err.println("ERROR: required single control file or separate file for each treatment.")
                    System.exit(1)
                }
                treatmentPaths.zip(controlPaths).map {
                    BinnedReadsQuery(genomeQuery, it.first, bin, it.second, unique = true, fragment = fragment)
                }
            }
        }
        return treatmentPaths.map { BinnedReadsQuery(genomeQuery, it, bin, unique = true, fragment = fragment) }
    }

    private fun configurePaths(outputPath: Path, chromSizesPath: Path) {
        if (ignoreConfigurePaths) {
            LOG.debug("IGNORE configurePaths")
            return
        }
        outputPath.createDirectories()
        Configuration.experimentsPath = outputPath
        Configuration.genomesPath = chromSizesPath.parent
        // See [Genome#chromSizesPath] for details
        System.getProperties().setProperty("chrom.sizes", chromSizesPath.toString())
    }

    private fun loadGenomeQuery(chromSizesPath: Path, chromosomes: List<String>): GenomeQuery {
        val fileName = chromSizesPath.fileName.toString()
        if (fileName.endsWith(".chrom.sizes")) {
            return GenomeQuery(fileName.substringBeforeLast(".chrom.sizes"), *chromosomes.toTypedArray())
        } else {
            throw IllegalArgumentException(
                    "ERROR: Unexpected reference name, please use the" +
                            "name of form <build>.chrom.sizes.")
        }
    }


    private fun getOptionParser(bin: Boolean = true,
                                fdr: Boolean = true,
                                gap: Boolean = true): OptionParser = object : OptionParser() {
        init {
            acceptsAll(listOf("d", "debug"), "Print all the debug information, used for troubleshooting.")
            acceptsAll(listOf("q", "quiet"), "Turn off output")
            acceptsAll(listOf("cs", "chrom.sizes"), "Chromosome sizes path, can be downloaded at\n" +
                    "http://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes")
                    .withRequiredArg().required()
                    .withValuesConvertedBy(PathConverter.exists())
            accepts("only").withOptionalArg().describedAs("chr1,chrN,...")
                    .withValuesSeparatedBy(',')

            acceptsAll(listOf("o", "output"), "Path to result bed")
                    .withRequiredArg()
                    .withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(listOf("fragment"), "Fragment size, read length if not given")
                    .withRequiredArg()
                    .ofType(Int::class.java)
            if (bin) {
                acceptsAll(listOf("b", "bin"), "Bin size")
                        .withRequiredArg()
                        .ofType(Int::class.java)
                        .defaultsTo(Span.BIN)
            }
            if (fdr) {
                acceptsAll(listOf("f", "fdr"), "Fdr value")
                        .withRequiredArg()
                        .ofType(Double::class.java)
                        .defaultsTo(Span.FDR)
            }
            if (gap) {
                acceptsAll(listOf("g", "gap"), "Gap size to merge peaks")
                        .withRequiredArg()
                        .ofType(Int::class.java)
                        .defaultsTo(Span.GAP)
            }
            acceptsAll(listOf("w", "workdir"), "Path to the working dir")
                    .withRequiredArg().withValuesConvertedBy(PathConverter.exists())
                    .defaultsTo(System.getProperty("user.dir").toPath())
            acceptsAll(listOf("threads"), "Parallelism level")
                    .withRequiredArg()
                    .ofType(Int::class.java)
            formatHelpWith(BuiltinHelpFormatter(200, 2))
        }
    }

    private fun configureLogging(quiet: Boolean, debug: Boolean) {
        if (quiet) {
            (LogManager.getCurrentLoggers().toList() + listOf(LogManager.getRootLogger())).forEach {
                (it as Logger).level = Level.OFF
            }

            val nullPrintStream = PrintStream(ByteStreams.nullOutputStream())
            System.setOut(nullPrintStream)
            System.setErr(nullPrintStream)
        }
        if (debug) {
            with(Logger.getRootLogger()) {
                level = Level.DEBUG
            }
        }
    }
}
