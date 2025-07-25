package org.jetbrains.bio.span

import joptsimple.OptionParser
import joptsimple.OptionSet
import org.jetbrains.bio.experiment.Configuration
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.coverage.AutoFragment
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.span.fit.AbstractSpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BIN
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FDR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_MAX_ITERATIONS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FIT_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.util.*
import org.jetbrains.bio.util.FileSize.Companion.GB
import org.slf4j.LoggerFactory
import java.nio.file.Path

/**
 * Tool for analyzing and comparing ChIP-Seq data.
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
     * This is a HACK.
     * Since [Configuration] allows configuring experimentsPath only once,
     * SpanCLA fails to set up correct working directory, if launched within the same process.
     */
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
analyze                         Peak calling
compare                         Differential peak calling
"""
    private const val ANALYZE = "analyze"
    private const val COMPARE = "compare"

    @JvmStatic
    fun main(args: Array<String>) {
        if (args.isEmpty()) {
            System.err.println("ERROR: No command given.")
            System.err.println(HELP)
        } else {
            when (args[0]) {
                ANALYZE -> SpanCLAAnalyze.analyze(args.copyOfRange(1, args.size))
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
        // [Shpynov] This is a hack: since -Xmx8G results in about 7G max memory available
        // The reason of such a discrepancy is the size of the garbage collector's survivor space.
        // https://stackoverflow.com/questions/23701207/why-do-xmx-and-runtime-maxmemory-not-agree
        if (maxMemory < 7 * GB) {
            LOG.warn(
                """
                    Recommended memory settings ${FileSize(8 * GB)} are not set. Current settings: ${FileSize(maxMemory)}.
                    Please use java memory option '-Xmx8G' to configure memory available for SPAN.""".trimIndent()
            )
        }
    }


    /**
     * Create [OptionParser] common for both [analyze] and [compare] procedures.
     */
    internal fun getOptionParser(): OptionParser = object : OptionParser() {
        init {
            acceptsAll(listOf("d", "debug"), "Print all the debug information, used for troubleshooting")
            acceptsAll(listOf("q", "quiet"), "Turn off output")
            acceptsAll(
                listOf("fmt", "format"),
                """Format of input file, supported formats: ${
                    ReadsFormat.values().joinToString(", ") { it.name }
                    }. 
                    Text format can be in zip or gzip archive""".trimIndent()
            ).withRequiredArg()
            acceptsAll(
                listOf("m", "model"),
                """
                    This option is used to specify SPAN model path.
                    Required for further semi-supervised peak calling.
                    If not provided, will be created in working directory
                    """.trimIndent()
            )
                .withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("l", "log"),
                "Path to log file. If not provided, will be created in working directory"
            )
                .withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())
            acceptsAll(
                listOf("cs", "chrom.sizes"),
                """
                    Chromosome sizes path, can be downloaded at
                    https://hgdownload.cse.ucsc.edu/goldenPath/<build>/bigZips/<build>.chrom.sizes
                    """.trimIndent()
            ).requiredUnless("model").withRequiredArg().withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("chr", "chromosomes"),
                "Chromosome names to process, should be a comma-separated list"
            ).withRequiredArg()

            accepts(
                "clip",
                "Clip max threshold for fine-tune boundaries of according to the local signal, 0 to disable"
            )
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_CLIP_MAX_SIGNAL)

            acceptsAll(
                listOf("p", "peaks"),
                """
                    Resulting peaks file in ENCODE broadPeak (BED 6+3) format.
                    If omitted, only the model fitting step is performed
                    """.trimIndent()
            ).withRequiredArg().withValuesConvertedBy(PathConverter.noCheck())

            acceptsAll(
                listOf("fragment"),
                """
                    Fragment size. If provided, reads are shifted appropriately.
                    --fragment 0 recommended for ATAC-Seq data processing.
                    If not provided or auto, the shift is estimated from the data (default: auto)
                    """.trimIndent()
            ).withRequiredArg().withValuesConvertedBy(FragmentConverter())
            acceptsAll(listOf("b", "bin"), "Bin size (default: $SPAN_DEFAULT_BIN)")
                .withRequiredArg()
                .ofType(Int::class.java)

            acceptsAll(listOf("blacklist", "bs"), "Blacklist regions file")
                .withRequiredArg()
                .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(
                listOf("multiple"),
                "Multiple hypothesis testing adjustment"
            )
                .availableIf("peaks")
                .withRequiredArg()
                .ofType(String::class.java)
                .defaultsTo(SPAN_DEFAULT_MULTIPLE_TEST_CORRECTION.name)

            acceptsAll(
                listOf("f", "fdr"),
                "False Discovery Rate cutoff to call significant regions"
            )
                .availableIf("peaks")
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_FDR)

            acceptsAll(
                listOf("sensitivity"),
                "Configures background sensitivity for peaks. Estimated from the data"
            )
                .availableIf("peaks")
                .withRequiredArg()
                .ofType(Double::class.java)
            accepts("summits", "Report summits, recommended for ATAC-seq or single-cell ATAC-seq experiments")
            accepts(
                "gap",
                "Gap is used for broad histone marks to join adjacent peaks. Estimated from the data"
            )
                .availableIf("peaks")
                .availableUnless("summits")
                .withRequiredArg()
                .ofType(Int::class.java)


            acceptsAll(
                listOf("w", "workdir"),
                "Path to the working directory. Used to save coverage and model cache"
            )
                .withRequiredArg().withValuesConvertedBy(PathConverter.exists())
                .defaultsTo(System.getProperty("user.dir").toPath())
            acceptsAll(
                listOf("threads"),
                "Parallelism level (default: ${Runtime.getRuntime().availableProcessors()})"
            )
                .withRequiredArg()
                .ofType(Int::class.java)
            acceptsAll(
                listOf("i", "iterations"),
                "Maximum number of iterations for Expectation Maximisation (EM) algorithm"
            )
                .withRequiredArg()
                .ofType(Int::class.java)
                .defaultsTo(SPAN_DEFAULT_FIT_MAX_ITERATIONS)
            acceptsAll(
                listOf("tr", "threshold"),
                "Convergence threshold for EM algorithm, use --debug option to see detailed info"
            )
                .withRequiredArg()
                .ofType(Double::class.java)
                .defaultsTo(SPAN_DEFAULT_FIT_THRESHOLD)
            acceptsAll(
                listOf("kd", "keep-duplicates"),
                """
                    Keep duplicates.
                    By default, SPAN filters out redundant reads aligned at the same genomic position.
                    Recommended for bulk single cell ATAC-Seq data processing
                    """.trimIndent()
            )
            acceptsAll(
                listOf("kc", "keep-cache"),
                """
                    Keep cache files.
                    By default SPAN creates cache files in working directory and removes them after computation is done
                    """.trimIndent()
            )
        }
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
            "$propertyName ($fitInfoValue) differs from the command line argument ($commandLineValue)"
        }
        val property = commandLineValue ?: fitInfoValue ?: default
        if (log) {
            LOG.info("$propertyId: $property")
        }
        return property
    }

    internal fun readsFormat(options: OptionSet, log: Boolean): ReadsFormat? {
        var explicitFormat: ReadsFormat? = null
        if ("format" in options) {
            explicitFormat = ReadsFormat.valueOf(options.valueOf("format") as String)
            if (log) LOG.info("FORMAT: $explicitFormat")
        } else {
            if (log) LOG.info("FORMAT: auto")
        }
        return explicitFormat
    }

    internal fun getUnique(
        options: OptionSet, fitInformation: AbstractSpanAnalyzeFitInformation? = null, log: Boolean = false
    ) = !getProperty(
        "keep-duplicates" in options, fitInformation?.unique?.not(), false,
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

    internal fun getFitThreshold(
        options: OptionSet, log: Boolean = false
    ) = getProperty(
        options.valueOf("threshold") as Double?, null, SPAN_DEFAULT_FIT_THRESHOLD,
        "convergence threshold", "CONVERGENCE THRESHOLD", log
    )

    internal fun getFitMaxIteration(
        options: OptionSet, log: Boolean = false
    ) = getProperty(
        options.valueOf("iterations") as Int?, null, SPAN_DEFAULT_FIT_MAX_ITERATIONS,
        "max iterations", "MAX ITERATIONS", log
    )

    /**
     * Checks if the genome build stored in the fitInformation matches the build inferred from the chrom.sizes file.
     * Also checks if all the chromosomes are present in the chrom.sizes file and have the same length.
     */
    internal fun checkGenomeInFitInformation(
        chromSizesPath: Path,
        fitInformation: AbstractSpanAnalyzeFitInformation
    ) {
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

    internal fun matchTreatmentsAndControls(
        commandLineTreatmentPaths: List<Path>,
        commandLineControlPaths: List<Path>
    ): List<SpanDataPaths>? {
        if (commandLineTreatmentPaths.isEmpty()) {
            return null
        }
        val spreadControlPaths = if (commandLineControlPaths.isNotEmpty()) {
            when (commandLineControlPaths.size) {
                1 -> Array(commandLineTreatmentPaths.size) { commandLineControlPaths.single() }.toList()
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

