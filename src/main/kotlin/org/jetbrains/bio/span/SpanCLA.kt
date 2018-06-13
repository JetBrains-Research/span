package org.jetbrains.bio.span

import com.google.common.io.ByteStreams
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
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.util.*
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
    private val LOG = Logger.getLogger(SpanCLA::class.java)

    // Load build properties
    init {
        val resource = SpanCLA::class.java.getResource("/span-build.properties")
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
compare                         Differential peak calling mode (experimental, use at your own risk)
"""
    private const val ANALYZE = "analyze"
    private const val COMPARE = "compare"

    @JvmStatic
    fun main(args: Array<String>) {
        Logs.addConsoleAppender(Level.INFO)
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
                    "ChIP-seq treatment file. bam, bed, or .bed.gz file;\n" +
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
                configureLogging("quiet" in options, "debug" in options)
                configureParallelism(options.valueOf("threads") as Int?)

                val workingDir = options.valueOf("workdir") as Path
                println("WORKING DIR: $workingDir")
                println("THREADS: ${parallelismLevel()}")

                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                val chromosomes = options.valuesOf("only") as List<String>
                val bin = options.valueOf("bin") as Int
                val fragment = options.valueOf("fragment") as Int?
                val gap = options.valueOf("gap") as Int

                // We would like to reuse as much caching as possible, resolve all the symbolic links
                val treatmentPaths = (options.valuesOf("treatment") as List<Path>).map { it.toRealPath() }
                val controlPaths = (options.valuesOf("control") as List<Path>?)?.map { it.toRealPath() }

                val genomeQuery = loadGenomeQuery(chromSizesPath, chromosomes)

                val labelsPath = options.valueOf("labels") as Path?
                val outputBed = options.valueOf("output") as Path?

                val fdr = options.valueOf("fdr") as Double
                println("TREATMENT: ${treatmentPaths.joinToString(", ", transform = Path::toString)}")
                if (controlPaths != null && controlPaths.isNotEmpty()) {
                    println("CONTROL: ${controlPaths.joinToString(", ", transform = Path::toString)}")
                } else {
                    println("CONTROL: none")
                }
                println("CHROM.SIZES: $chromSizesPath")
                println("GENOME: ${genomeQuery.id}")
                println("FRAGMENT: $fragment")
                if (outputBed != null) {
                    if (labelsPath != null) {
                        println("LABELS: $labelsPath")
                        println("BIN, FDR, GAP options are ignored.")
                    } else {
                        println("BIN: $bin")
                        println("FDR: $fdr")
                        println("GAP: $gap")
                    }
                    println("OUTPUT: $outputBed")
                } else {
                    println("BIN: $bin")
                    println("NO output path given, process model fitting only.")
                    println("LABELS, FDR, GAP options are ignored.")
                }

                configurePaths(workingDir, chromSizesPath)
                val coverageQueries = coverageQueries(genomeQuery, treatmentPaths, controlPaths, fragment)

                val peakCallingExperiment = Span.getPeakCallingExperiment(genomeQuery, coverageQueries, bin)
                if (outputBed != null) {
                    if (labelsPath == null) {
                        val peaks = peakCallingExperiment.getPeaks(fdr, gap)
                        savePeaks(peaks, outputBed, "peak${if (fragment != null) "_$fragment" else ""}_${bin}_${fdr}_${gap}")
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
                        results.saveOptimalResults(outputBed.parent / "${outputBed.fileName.stem}_parameters.csv")
                        val peaks = Span.getPeakCallingExperiment(genomeQuery,
                                coverageQueries, SPAN.DEFAULT_BIN).getPeaks(optimalFDR, optimalGap)
                        savePeaks(peaks, outputBed, "peak${if (fragment != null) "_$fragment" else ""}_" +
                                "${SPAN.DEFAULT_BIN}_${optimalFDR}_$optimalGap")
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


    private fun compare(params: Array<String>) {
        with(getOptionParser()) {
            acceptsAll(listOf("t1", "treatment1"),
                    "ChIP-seq treatment file 1. bam, bed, or .bed.gz file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c1", "control1"),
                    "Control file 1. bam, bed or bed.gz file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(listOf("t2", "treatment2"),
                    "ChIP-seq treatment file 2. bam, bed, or .bed.gz file;\n" +
                            "If multiple files are given, treated as replicates.")
                    .withRequiredArg().required()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(listOf("c2", "control2"),
                    "Control file 2. bam, bed or bed.gz file;\n" +
                            "Single control file or separate file per each\n" +
                            "treatment file required.")
                    .withRequiredArg()
                    .withValuesSeparatedBy(",")
                    .withValuesConvertedBy(PathConverter.exists())

            parse(params) { options ->
                configureLogging("quiet" in options, "debug" in options)
                configureParallelism(options.valueOf("threads") as Int?)

                val workingDir = options.valueOf("workdir") as Path
                println("WORKING DIR: $workingDir")
                println("THREADS: ${parallelismLevel()}")

                val chromSizesPath = options.valueOf("chrom.sizes") as Path
                val chromosomes = options.valuesOf("only") as List<String>

                val fragment = options.valueOf("fragment") as Int?
                val bin = options.valueOf("bin") as Int
                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double

                val treatmentPaths1 = options.valuesOf("treatment1") as List<Path>
                val treatmentPaths2 = options.valuesOf("treatment2") as List<Path>
                val controlPaths1 = options.valuesOf("control1") as List<Path>?
                val controlPaths2 = options.valuesOf("control2") as List<Path>?
                val outputBed = options.valueOf("output") as Path?
                val genomeQuery = loadGenomeQuery(chromSizesPath, chromosomes)

                println("TREATMENT1: ${treatmentPaths1.joinToString(", ", transform = Path::toString)}")
                if (controlPaths1 != null && controlPaths1.isNotEmpty()) {
                    println("CONTROL1: ${controlPaths1.joinToString(", ", transform = Path::toString)}")
                } else {
                    println("CONTROL1: none")
                }
                println("TREATMENT2: ${treatmentPaths2.joinToString(", ", transform = Path::toString)}")
                if (controlPaths2 != null && controlPaths2.isNotEmpty()) {
                    println("CONTROL2: ${controlPaths2.joinToString(", ", transform = Path::toString)}")
                } else {
                    println("CONTROL2: none")
                }
                println("CHROM.SIZES: $chromSizesPath")
                println("GENOME: ${genomeQuery.id}")
                println("FRAGMENT: $fragment")
                println("BIN: $bin")
                println("FDR: $fdr")
                println("GAP: $gap")
                if (outputBed != null) {
                    println("OUTPUT: $outputBed")
                } else {
                    println("BIN: $bin")
                    println("NO output path given, process model fitting only.")
                    println("LABELS, FDR, GAP options are ignored.")
                }
                configurePaths(workingDir, chromSizesPath)

                val coverageQueries1 = coverageQueries(genomeQuery, treatmentPaths1, controlPaths1, fragment)
                val coverageQueries2 = coverageQueries(genomeQuery, treatmentPaths2, controlPaths2, fragment)
                val experiment = Span.getDifferentialPeakCallingExperiment(
                        genomeQuery, coverageQueries1, coverageQueries2, bin)
                if (outputBed != null) {
                    val peaks = experiment.getPeaks(fdr, gap)
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
                                fragment: Int?): List<ReadsQuery> {
        if (controlPaths != null && controlPaths.isNotEmpty()) {
            if (controlPaths.size == 1) {
                return treatmentPaths.map {
                    ReadsQuery(genomeQuery, it, controlPaths.first(), unique = true, fragment = fragment)
                }
            } else {
                if (controlPaths.size != treatmentPaths.size) {
                    System.err.println("ERROR: required single control file or separate file for each treatment.")
                    System.exit(1)
                }
                treatmentPaths.zip(controlPaths).map {
                    ReadsQuery(genomeQuery, it.first, it.second, unique = true, fragment = fragment)
                }
            }
        }
        return treatmentPaths.map { ReadsQuery(genomeQuery, it, unique = true, fragment = fragment) }
    }

    private fun configurePaths(outputPath: Path, chromSizesPath: Path) {
        outputPath.createDirectories()
        Configuration.setWorkDir(outputPath)
        System.getProperties().setProperty("chrom.sizes", chromSizesPath.toString())
        System.getProperties().setProperty("genomes.path", chromSizesPath.parent.toString())
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
            acceptsAll(listOf("cs", "chrom.sizes"), "Chromosome sizes path, can be downloaded at " +
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
