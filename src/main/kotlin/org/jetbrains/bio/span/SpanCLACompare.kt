package org.jetbrains.bio.span

import joptsimple.OptionSet
import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanDifferentialPeakCallingExperiment
import org.jetbrains.bio.experiments.fit.SpanFitResults
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.peaks.getFdrGapPeaks
import org.jetbrains.bio.util.*
import org.slf4j.event.Level
import java.nio.file.Path

object SpanCLACompare {

    internal fun compare(params: Array<String>) {
        with(SpanCLA.getOptionParser()) {
            acceptsAll(
                listOf("t1", "treatment1"),
                "ChIP-seq treatment file 1. bam, bed or .bed.gz file;\n" +
                        "If multiple files are given, treated as replicates."
            )
                .withRequiredArg().required()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("c1", "control1"),
                "Control file 1. bam, bed or .bed.gz file;\n" +
                        "Single control file or separate file per each\n" +
                        "treatment file required."
            )
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())

            acceptsAll(
                listOf("t2", "treatment2"),
                "ChIP-seq treatment file 2. bam, bed or .bed.gz file;\n" +
                        "If multiple files are given, treated as replicates."
            )
                .withRequiredArg().required()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())
            acceptsAll(
                listOf("c2", "control2"),
                "Control file 2. bam, bed or .bed.gz file;\n" +
                        "Single control file or separate file per each\n" +
                        "treatment file required."
            )
                .withRequiredArg()
                .withValuesSeparatedBy(",")
                .withValuesConvertedBy(PathConverter.exists())


            parse(params) { options ->
                if ("quiet" in options) {
                    Logs.quiet()
                } else {
                    Logs.addConsoleAppender(if ("debug" in options) Level.DEBUG else Level.INFO)
                }
                SpanCLA.LOG.info("SPAN ${SpanCLA.version()}")
                SpanCLA.LOG.info("COMMAND: compare ${params.joinToString(" ")}")

                val workingDir = options.valueOf("workdir") as Path

                val gap = options.valueOf("gap") as Int
                val fdr = options.valueOf("fdr") as Double
                val peaksType = options.valueOf("peaks-type") as String
                val peaksPath = options.valueOf("peaks") as Path?
                val threads = options.valueOf("threads") as Int?

                // Configure logging
                val fragment = SpanCLA.getFragment(options)
                val binSize = SpanCLA.getBin(options)
                val id = if (peaksPath != null) {
                    peaksPath.stemGz
                } else {
                    // No peaks, generate ID from command-line options.
                    // Option parser guarantees that treatment paths are not empty here.
                    val datas = getComparePaths(options)
                    val ids = datas.toList().flatMap { data ->
                        listOf(data.map { it.treatment }, data.mapNotNull { it.control }).flatMap { paths ->
                            paths.map { it.stemGz }
                        }
                    }.toMutableList()
                    ids.add(binSize.toString())
                    if (fragment is FixedFragment) {
                        ids.add(fragment.size.toString())
                    }
                    val unique = SpanCLA.getUnique(options)
                    if (unique) {
                        ids.add("unique")
                    }
                    reduceIds(ids)
                }

                val logPath = SpanCLA.configureLogFile(workingDir, id)
                SpanCLA.LOG.info("LOG: $logPath")

                // Call now to preserve correct params logging
                val lazyDifferentialPeakCallingResults = differentialPeakCallingResults(options)

                SpanCLA.LOG.info("FDR: $fdr")
                SpanCLA.LOG.info("GAP: $gap")
                if (peaksPath != null) {
                    SpanCLA.LOG.info("TYPE: $peaksType")
                    SpanCLA.LOG.info("PEAKS: $peaksPath")
                } else {
                    SpanCLA.LOG.info("NO output path given, process model fitting only.")
                    SpanCLA.LOG.info("LABELS, FDR, GAP options are ignored.")
                }
                configureParallelism(threads)
                SpanCLA.LOG.info("THREADS: ${parallelismLevel()}")

                SpanCLA.checkMemory()

                val differentialPeakCallingResults = lazyDifferentialPeakCallingResults.value
                val genomeQuery = differentialPeakCallingResults.fitInfo.genomeQuery()
                if (peaksPath != null) {
                    val peaks = differentialPeakCallingResults.getFdrGapPeaks(genomeQuery, fdr, gap)
                    Peak.savePeaks(
                        peaks, peaksPath,
                        "diff${if (fragment is FixedFragment) "_$fragment" else ""}_${binSize}_${fdr}_${gap}"
                    )
                    SpanCLA.LOG.info("Saved result to $peaksPath")
                }
            }
        }
    }


    /**
     * Retrieves the paths (treatment1, optional control1), (treatment2, optional control2)
     * Checks that they are consistent.
     */
    private fun getComparePaths(
        options: OptionSet,
        log: Boolean = false
    ): Pair<List<SpanDataPaths>, List<SpanDataPaths>> {
        val treatmentPaths1 = options.valuesOf("treatment1") as List<Path>
        val treatmentPaths2 = options.valuesOf("treatment2") as List<Path>
        val controlPaths1 = options.valuesOf("control1") as List<Path>
        val controlPaths2 = options.valuesOf("control2") as List<Path>

        val paths1 = SpanCLA.getCommandLinePaths(treatmentPaths1, controlPaths1)
        check(paths1 != null) { "No treatment files provided for set 1, use -t1 option." }
        val paths2 = SpanCLA.getCommandLinePaths(treatmentPaths2, controlPaths2)
        check(paths2 != null) { "No treatment files provided for set 2, use -t2 option." }

        if (log) {
            SpanCLA.LOG.info("TREATMENT1: ${treatmentPaths1.joinToString(", ", transform = Path::toString)}")
            if (controlPaths1.isNotEmpty()) {
                SpanCLA.LOG.info("CONTROL1: ${controlPaths1.joinToString(", ", transform = Path::toString)}")
            } else {
                SpanCLA.LOG.info("CONTROL1: none")
            }
            SpanCLA.LOG.info("TREATMENT2: ${treatmentPaths2.joinToString(", ", transform = Path::toString)}")
            if (controlPaths2.isNotEmpty()) {
                SpanCLA.LOG.info("CONTROL2: ${controlPaths2.joinToString(", ", transform = Path::toString)}")
            } else {
                SpanCLA.LOG.info("CONTROL2: none")
            }
        }
        return paths1 to paths2
    }

    /**
     * Configure logging and get [SpanFitResults] in a most concise and effective way.
     * Parses and logs most of the command line arguments.
     */
    private fun differentialPeakCallingResults(options: OptionSet): Lazy<SpanFitResults> {
        val chromSizesPath = SpanCLA.getAndLogWorkDirAndChromSizes(options)
        val genomeQuery = GenomeQuery(Genome[chromSizesPath!!])
        val (data1, data2) = getComparePaths(options, log = true)
        SpanCLA.LOG.info("CHROM.SIZES: $chromSizesPath")
        val bin = SpanCLA.getBin(options, log = true)
        val fragment = SpanCLA.getFragment(options, log = true)
        val unique = SpanCLA.getUnique(options, log = true)
        val maxIter = SpanCLA.getMaxIter(options, log = true)
        val threshold = SpanCLA.getThreshold(options, log = true)
        return lazy {
            val experiment = SpanDifferentialPeakCallingExperiment.getExperiment(
                genomeQuery, data1, data2, bin, fragment, unique,
                threshold, maxIter
            )
            experiment.results
        }
    }

}