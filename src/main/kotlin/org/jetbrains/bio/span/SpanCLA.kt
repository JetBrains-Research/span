package org.jetbrains.bio.span

import com.google.common.annotations.VisibleForTesting
import joptsimple.OptionParser
import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.fit.AbstractSpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_BIN
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_FDR
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_GAP
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.util.*
import org.jetbrains.bio.util.FileSize.Companion.GB
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * Tool for analyzing and comparing ChIP-Seq data.
 * `analyze` and `compare` procedures rely on the Zero Inflated Negative Binomial Restricted Algorithm.
 *
 * See [SpanCLAAnalyze] for peak calling procedure.
 * See [SpanCLACompare] for differential peak calling.
 *
 * @author Oleg Shpynov
 * @since  14/09/15
 */
@Suppress("UNCHECKED_CAST")
object SpanCLA {
    internal val LOG = LoggerFactory.getLogger(SpanCLA::class.java)

    /**
     * Shpynov:
     * Since [Configuration] allows configuring experimentsPath only once,
     * SpanCLA fails to set up correct working directory, if launched within the same process.
     * This is a HACK.
     */
    @VisibleForTesting
    var ignoreConfigurePaths: Boolean = false

    init {
        // Load build properties
        val resource = SpanCLA::class.java.getResource("/span.properties")
        resource?.openStream()?.use { System.getProperties().load(it) }
    }

    internal fun version() =
        "${System.getProperty("span.build.version", "@VERSION@.@build@")} " +
                "built on ${System.getProperty("span.build.date", "@DATE@")}"


    private const val HELP = """
Option                          Description
---------------------           -----------
-?, -h, --help                  Show help
-v, --version                   Show version
analyze                         Peak calling mode
compare                         Differential peak calling mode 
experimental                    Experimental features
"""
    private const val ANALYZE = "analyze"
    private const val EXPERIMENTAL = "experimental"
    private const val COMPARE = "compare"

    @JvmStatic
    fun main(args: Array<String>) {
        if (args.isEmpty()) {
            System.err.println("ERROR: No command given.")
            System.err.println(HELP)
        } else {
            when (args[0]) {
                ANALYZE -> SpanCLAAnalyze.analyze(args.copyOfRange(1, args.size), false)
                EXPERIMENTAL -> SpanCLAAnalyze.analyze(args.copyOfRange(1, args.size), true)
                COMPARE -> SpanCLACompare.compare(args.copyOfRange(1, args.size))

                "-?", "-h", "--help" -> println(HELP)
                "-v", "--version" -> println(version())

                else -> {
                    System.err.println("ERROR: Unknown command: ${args[0]}.")
                    System.err.println(HELP)
                }
            }
        }
    }


    internal fun checkMemory() {
        val maxMemory = Runtime.getRuntime().maxMemory()
        // [Shpynov] This is a hack: since -Xmx4G results in about 3.5G max memory available
        // The reason of such a discrepancy is the size of the garbage collector's survivor space.
        // https://stackoverflow.com/questions/23701207/why-do-xmx-and-runtime-maxmemory-not-agree
        // 3.5 works, however use decreased level to be sure.
        if (maxMemory < 3 * GB) {
            LOG.warn(
                "Recommended memory settings ${FileSize(4 * GB)} are not set. Current settings: ${FileSize(maxMemory)}.\n" +
                        "Please use java memory option '-Xmx4G' to configure memory available for SPAN."
            )
        }
    }


    /**
     * Initialize [Configuration] with given paths
     */
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


    /**
     * Create [OptionParser] common for both [analyze] and [compare] procedures.
     */
    internal fun getOptionParser(): OptionParser = object : OptionParser() {
        init {
            acceptsAll(listOf("d", "debug"), "Print all the debug information, used for troubleshooting.")
            acceptsAll(listOf("q", "quiet"), "Turn off output")
            acceptsAll(listOf("m", "model"), "Path to model file")
                .withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("cs", "chrom.sizes"),
                "Chromosome sizes path, can be downloaded at\n" +
                        "https://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes"
            ).requiredUnless("model").withRequiredArg().withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("p", "peaks"), "Path to result peaks file in ENCODE broadPeak (BED 6+3) format"
            ).withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("fragment"),
                "Fragment size. If it's an integer, reads are shifted appropriately.\n" +
                        "If it's the string 'auto', the shift is estimated from the data. (default: auto)"
            ).withRequiredArg().withValuesConvertedBy(FragmentConverter())
            acceptsAll(listOf("b", "bin"), "Bin size. (default: ${SPAN_DEFAULT_BIN})")
                .withRequiredArg()
                .ofType(Int::class.java)

            acceptsAll(listOf("f", "fdr"), "FDR value.")
                .availableIf("peaks")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_FDR)
            acceptsAll(listOf("g", "gap"), "Gap size to merge consequent peaks.")
                .availableIf("peaks")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(SPAN_DEFAULT_GAP)
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
            acceptsAll(listOf("i", "iterations"), "Maximum number of iterations for EM algorithm")
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(Fitter.MAX_ITERATIONS)
            acceptsAll(
                listOf("tr", "threshold"), "Convergence threshold for EM algorithm, " +
                        "use --debug option to see detailed info"
            )
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(Fitter.THRESHOLD)
        }
    }

    /**
     * Configure logging to file. Returns path to log file.
     */
    internal fun configureLogFile(workingDir: Path, id: String): Path =
        (workingDir / "logs" / "$id.log").apply {
            Logs.addLoggingToFile(this)
        }

    /**
     * Compare the command-line option (null if not given) and stored value (null if not present).
     * Fail if these differ. Return default if both are null. Log the value if [log] is true.
     */
    internal fun <T> getProperty(
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

    internal fun getUnique(
        options: OptionSet, fitInformation: AbstractSpanAnalyzeFitInformation? = null, log: Boolean = false
    ) = !getProperty(
        if ("keep-dup" in options) options.valueOf("keep-dup") as Boolean else null,
        fitInformation?.unique?.not(), false,
        "'keep duplicates' flag", "KEEP DUPLICATES", log
    )

    internal fun getFragment(
        options: OptionSet, fitInformation: AbstractSpanAnalyzeFitInformation? = null, log: Boolean = false
    ) = getProperty(
        options.valueOf("fragment") as Fragment?, fitInformation?.fragment, AutoFragment,
        "fragment size", "FRAGMENT", log
    )

    internal fun getBin(
        options: OptionSet, fitInformation: SpanFitInformation? = null, log: Boolean = false
    ) = getProperty(
        options.valueOf("bin") as Int?, fitInformation?.binSize, SPAN_DEFAULT_BIN,
        "bin size", "BIN", log
    )

    internal fun getThreshold(
        options: OptionSet, log: Boolean = false
    ) = getProperty(
        options.valueOf("threshold") as Double?, null, Fitter.THRESHOLD,
        "convergence threshold", "CONVERGENCE THRESHOLD", log
    )

    internal fun getMaxIter(
        options: OptionSet, log: Boolean = false
    ) = getProperty(
        options.valueOf("iterations") as Int?, null, Fitter.MAX_ITERATIONS,
        "max iterations", "MAX ITERATIONS", log
    )

    internal fun getAndLogWorkDirAndChromSizes(
        options: OptionSet, fitInformation: AbstractSpanAnalyzeFitInformation? = null
    ): Path? {
        val workingDir = options.valueOf("workdir") as Path
        LOG.info("WORKING DIR: $workingDir")

        val chromSizesPath = options.valueOf("chrom.sizes") as Path?
        configurePaths(workingDir, chromSizesPath)

        if (fitInformation != null && chromSizesPath != null) {
            val genome = Genome[chromSizesPath]
            check(genome.build == fitInformation.build) {
                "Stored genome build ${genome.build} differs from " +
                        "the one inferred from the chrom.sizes file $chromSizesPath"
            }
            val chromosomeMap = genome.chromosomeNamesMap
            // we don't check the map equality, since the stored map contains only non-empty chromosomes
            fitInformation.chromosomesSizes.forEach { (name, length) ->
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
        return chromSizesPath
    }

    internal fun getCommandLinePaths(
        commandLineTreatmentPaths: List<Path>,
        commandLineControlPaths: List<Path>
    ): List<SpanDataPaths>? {
        if (commandLineTreatmentPaths.isEmpty()) {
            return null
        }
        val spreadControlPaths = if (commandLineControlPaths.isNotEmpty()) {
            when (commandLineControlPaths.size) {
                1 ->
                    Array(commandLineTreatmentPaths.size) { commandLineControlPaths.single() }.toList()
                commandLineTreatmentPaths.size -> commandLineControlPaths
                else -> throw IllegalArgumentException(
                    "Expected either: no control files, a single control file, " +
                            "or as many control files as treatment files."
                )
            }
        } else {
            Array(commandLineTreatmentPaths.size) { null as Path? }.toList()
        }
        return commandLineTreatmentPaths.zip(spreadControlPaths).map { (treatment, control) ->
            SpanDataPaths(treatment, control)
        }
    }

}

