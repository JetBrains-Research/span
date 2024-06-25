package org.jetbrains.bio.span

import com.google.common.collect.MinMaxPriorityQueue
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.unpack
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.AUTOCORRELATION_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_MODERATE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_SEVERE
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.ModelToPeaks.computeCorrelations
import org.jetbrains.bio.span.peaks.ModelToPeaks.detectGapModel
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateCandidatesNumberLenDist
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateMinPivotGap
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateSensitivityGap
import org.jetbrains.bio.span.peaks.ModelToPeaks.getChromosomeCandidates
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_BACKGROUND_SENSITIVITY_VARIANTS
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_GAPS_VARIANTS
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.io.BufferedWriter
import java.io.FileWriter
import java.nio.file.Path
import java.util.*
import kotlin.math.*

object SpanResultsAnalysis {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun doDeepAnalysis(
        actualModelPath: Path,
        spanFitResults: SpanFitResults,
        fitInfo: SpanAnalyzeFitInformation,
        genomeQuery: GenomeQuery,
        fdr: Double,
        blacklistPath: Path?,
        peaksList: List<Peak>,
        peaksPath: Path?
    ) {
        check(fitInfo.normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        check(fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }) {
            "Coverage information is not available"
        }
        val infoFile = if (peaksPath != null) "$peaksPath.txt" else null
        infoFile?.toPath()?.deleteIfExists()
        val infoWriter = if (infoFile != null) BufferedWriter(FileWriter(infoFile)) else null

        // Save basic stats to infoFile
        LOG.info("Processing basic info")
        if (infoWriter != null) {
            val aboutModel = spanFitResults.modelInformation(actualModelPath)
            val aboutPeaks = PeaksInfo.compute(
                genomeQuery,
                peaksList.map { it.location }.stream(),
                peaksPath!!.toUri(),
                fitInfo.paths.map { it.treatment }
            )
            infoWriter.write((aboutModel + aboutPeaks).joinToString("\n") { (k, v) ->
                "${k.name}: ${k.render(v)}"
            } + "\n")
        }

        LOG.info("Analysing auto correlations...")
        val ncq = fitInfo.normalizedCoverageQueries!!.first()
        val (treatmentScale, controlScale, beta, minCorrelation) = ncq.coveragesNormalizedInfo
        val treatmentCoverage = ncq.treatmentReads.coverage()
        val treatmentTotal = genomeQuery.get().sumOf {
            treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange).toLong()
        }
        logInfo("Treatment coverage: $treatmentTotal", infoWriter)
        val controlCoverage = ncq.controlReads?.coverage()
        if (controlCoverage != null) {
            val controlTotal = genomeQuery.get().sumOf {
                controlCoverage.getBothStrandsCoverage(it.chromosomeRange).toLong()
            }
            logInfo("Treatment scale: $treatmentScale", infoWriter)
            logInfo("Control coverage: $controlTotal", infoWriter)
            logInfo("Control scale: $controlScale", infoWriter)
            logInfo("Beta: $beta", infoWriter)
            logInfo("Min control correlation: $minCorrelation", infoWriter)
        }

        val maxGapCoverage = detectGapCoverage(
            genomeQuery, treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta,
            fitInfo.binSize, blacklistPath
        )
        logInfo("Maximal gap coverage: $maxGapCoverage", infoWriter)

        val maxGapModel = detectGapModel(
            genomeQuery, spanFitResults, blacklistPath
        )
        logInfo("Maximal gap model: $maxGapModel", infoWriter)

        LOG.info("Analysing tracks roughness...")
        val roughness = detectAverageRoughness(
            genomeQuery, treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
        )
        logInfo("Track roughness: ${"%.3f".format(roughness)}", infoWriter)
        val roughnessNoControl = detectAverageRoughness(
            genomeQuery, treatmentCoverage, treatmentScale, null, null, 0.0
        )
        logInfo("Track roughness no control: ${"%.3f".format(roughnessNoControl)}", infoWriter)

        if (blacklistPath != null) {
            LOG.info("Computing blacklisted roughness")
            val roughnessBlackListed = detectAverageRoughness(
                genomeQuery, treatmentCoverage, 1.0, null, null, 0.0, blacklistPath
            )
            logInfo(
                "Track roughness blacklisted: ${"%.3f".format(roughnessBlackListed)}",
                infoWriter
            )
        }

        LOG.info("Analysing min pivot gap")
        val minPivotGap = estimateMinPivotGap(genomeQuery, spanFitResults, fdr)
        logInfo("Minimal pivot gap: $minPivotGap", infoWriter)
        if (minPivotGap != null) {
            when {
                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_SEVERE ->
                    logInfo("Quality: severe fragmentation", infoWriter)
                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_MODERATE ->
                    logInfo("Quality: moderate fragmentation", infoWriter)
                else ->
                    logInfo("Quality: light fragmentation", infoWriter)
            }
        } else {
            logInfo("Quality: good", infoWriter)
        }

        val (sensitivity2use, gap2use) = estimateSensitivityGap(
            genomeQuery,
            spanFitResults,
            null,
            null,
            minPivotGap
        )
        logInfo("Estimated sensitivity: $sensitivity2use", infoWriter)
        logInfo("Estimated gap: $gap2use", infoWriter)

        infoWriter?.close()

        LOG.info("Analysing candidates characteristics wrt sensitivity and gap...")
        val sensDetailsFile = if (peaksPath != null) "$peaksPath.sensitivity.tsv" else null
        sensDetailsFile?.toPath()?.deleteIfExists()
        val sensDetailsWriter = if (sensDetailsFile != null)
            BufferedWriter(FileWriter(sensDetailsFile))
        else
            null
        logInfo(
            "Sensitivity\tGap\tCandidatesN\tCandidatesAL\tCandidatesAD", sensDetailsWriter, false
        )
        for (bgs in SPAN_BACKGROUND_SENSITIVITY_VARIANTS) {
            for (g in SPAN_GAPS_VARIANTS) {
                val (candidatesN, candidatesAL, candidatesAD) =
                    estimateCandidatesNumberLenDist(genomeQuery, spanFitResults, fdr, bgs, g, null)
                logInfo(
                    "$bgs\t$g\t$candidatesN\t$candidatesAL\t$candidatesAD", sensDetailsWriter, false
                )
            }
        }
        sensDetailsWriter?.close()

        LOG.info("Analysing peaks segmentation wrt sensitivity")
        var bedtoolsPresent = true
        try {
            LOG.info("Checking if bedtools is installed...")
            Exec.exec("bedtools", "--help", output = OutputType.TEXT)
        } catch (e: Exception) {
            bedtoolsPresent = false
        }
        if (peaksPath != null) {
            if (!bedtoolsPresent) {
                LOG.warn("bedtools not available. Cannot create sensitivity track view.")
            } else {
                val blackList =
                    if (blacklistPath != null) LocationsMergingList.load(genomeQuery, blacklistPath) else null
                withTempDirectory(peaksPath.fileName!!.stem) { dir ->
                    LOG.info("Saving sensitivity peaks to $dir")
                    val peaksPaths = arrayListOf<String>()
                    val peaksSens = doubleArrayOf(10.0, 5.0, 1.0, 0.1, 1e-4, 1e-6)
                    for (s in peaksSens) {
                        val candidates = genomeQuery.get()
                            .flatMap { chr ->
                                if (!spanFitResults.fitInfo.containsChromosomeInfo(chr)) {
                                    return@flatMap emptyList<ChromosomeRange>()
                                }
                                getChromosomeCandidates(
                                    spanFitResults,
                                    chr,
                                    fdr,
                                    s,
                                    SPAN_DEFAULT_GAP
                                ).first.map { it.on(chr) }
                            }
                            .filter { blackList == null || !blackList.intersects(it.location) }
                        val path = dir / "peaks_$s.peak"
                        CSVFormat.TDF.print(path.bufferedWriter()).use { printer ->
                            candidates.forEach { r ->
                                printer.printRecord(r.chromosome.name, r.startOffset.toString(), r.endOffset.toString())
                            }
                        }
                        peaksPaths.add(path.toString())
                    }
                    val out = dir / "multiinter.peak"
                    val sensitivityTrack = "$peaksPath.sensitivity.bed"
                    sensitivityTrack.toPath().deleteIfExists()
                    withTempFile("script", ".sh") { script ->
                        script.bufferedWriter().use { w ->
                            w.write("bedtools multiinter -i ${peaksPaths.joinToString(" ")} | awk -v OFS='\\t' '{print $1, $2, $3, $5}' > $out")
                        }
                        "bash".exec(script.toString())
                    }
                    val scoredLocations = arrayListOf<XLocation<String>>()
                    val format = BedFormat.auto(out)
                    format.parse(out) {
                        it.forEach { entry ->
                            val chromosome = genomeQuery[entry.chrom]
                            if (chromosome != null) {
                                val e = entry.unpack(format)
                                scoredLocations.add(XLocation(Location(e.start, e.end, chromosome), e.name))
                            }
                        }
                    }
                    val lls = Array(peaksPaths.size) { LocationsMergingList.builder(genomeQuery) }
                    scoredLocations.sortedBy { it.l }.forEach { r ->
                        val idx = r.x.split(',').minOf { it.toInt() }
                        lls[idx - 1].add(r.l)
                    }
                    val peaks =
                        lls.flatMapIndexed { i, ll -> ll.build().toList().map { XLocation(it, i) } }.sortedBy { it.l }
                    CSVFormat.TDF.print(sensitivityTrack.toPath().bufferedWriter()).use { printer ->
                        peaks.forEach { r ->
                            val score = (1000.0 * (peaksSens.size - r.x) / peaksSens.size).toInt()
                            val s = peaksSens[r.x]
                            printer.printRecord(
                                r.l.chromosome.name,
                                r.l.startOffset * fitInfo.binSize,
                                r.l.endOffset * fitInfo.binSize,
                                "s$s",
                                score
                            )
                        }
                        LOG.info("Sensitivity file saved to $sensitivityTrack")
                    }
                }
            }
        }
    }

    private fun logInfo(msg: String, infoWriter: BufferedWriter?, useLog: Boolean = true) {
        if (useLog)
            SpanCLA.LOG.info(msg)
        else
            println(msg)
        if (infoWriter != null) {
            infoWriter.write(msg)
            infoWriter.newLine()
        }
    }

    class XLocation<T>(val l: Location, val x: T) : LocationAware {
        override val location: Location
            get() = l
    }

    val REGION_LEN = 10_000
    val TOP_REGIONS = 10_000
    val WORK_REGIONS = 200
    val RESOLUTION = 100

    /**
     * Detects roughness of the coverage.
     */
    fun detectAverageRoughness(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        treatmentScale: Double,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        blackListPath: Path? = null,
        region_len: Int = REGION_LEN,
        top_regions: Int = TOP_REGIONS,
        work_regions: Int = WORK_REGIONS,
        resolution: Int = RESOLUTION
    ): Double {
        LOG.debug("Compute coverage in regions")
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange) > 0 }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val genome_regions = limitedQuery.get().sumOf { floor(it.length.toDouble() / region_len).toLong() }
        check(top_regions < genome_regions) {
            "Too many top regions $top_regions > $genome_regions"
        }

        val comparator = Comparator<Triple<Chromosome, Int, Double>> { o1, o2 -> o2.third.compareTo(o1.third) }
        val region_coverages: MinMaxPriorityQueue<Triple<Chromosome, Int, Double>> =
            MinMaxPriorityQueue
                .orderedBy(comparator)
                .maximumSize(top_regions)
                .create()
        var regions = 0
        val blackList = if (blackListPath != null) LocationsMergingList.load(limitedQuery, blackListPath) else null
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            for (i in 0 until floor(chr.length.toDouble() / region_len).toInt()) {
                regions += 1
                val start = region_len * i
                val end = region_len * (i + 1)
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                val c = coverage(
                    chr, start, end,
                    treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
                )
                region_coverages.add(Triple(chr, i, c))
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, regions)
        }

        val region_coverages_array = region_coverages.toTypedArray()
        region_coverages_array.sortWith(comparator)

        val step = if (region_coverages_array.size > work_regions) {
            LOG.debug("Pick $work_regions / $top_regions uniform regions for computation speedup")
            ceil(region_coverages_array.size.toDouble() / work_regions).toInt()
        } else
            1

        val std_means = DoubleArray(work_regions)
        for (i in region_coverages_array.indices) {
            if (i % step != 0) {
                continue
            }
            val (chr, start, _) = region_coverages_array[i]
            val stats = DoubleArray(region_len / resolution) {
                coverage(
                    chr, start + it * resolution, start + (it + 1) * resolution,
                    treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
                )
            }
            val mean = stats.average()
            val std = stats.standardDeviation()
            std_means[i / step] = if (mean > 0) std / mean else 0.0
        }
        return std_means.average()
    }


    /**
     * Detects maximal gap, that signal autocorrelation is >= minCorrelation
     */
    fun detectGapCoverage(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        treatmentScale: Double,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        bin: Int,
        blackListPath: Path? = null,
        minCorrelation: Double = AUTOCORRELATION_THRESHOLD
    ): Int {
        LOG.debug("Compute maximal gap by coverage")
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange) > 0 }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val totalBins = limitedQuery.get().sumOf { floor(it.length.toDouble() / bin).toInt() }
        val coverage = DoubleArray(totalBins) { 0.0 }
        var i = 0
        val blackList = if (blackListPath != null) LocationsMergingList.load(limitedQuery, blackListPath) else null
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            for (j in 0 until floor(chr.length.toDouble() / bin).toInt()) {
                val start = j * bin
                val end = (j + 1) * bin
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                coverage[i++] = coverage(
                    chr, start, end,
                    treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
                )
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, totalBins)
        }
        val correlations = computeCorrelations(coverage)
        val maxGap = correlations.indices.filter { correlations[it] >= minCorrelation }.maxOrNull()!!
        return maxGap
    }


    private fun coverage(
        chromosome: Chromosome,
        start: Int, end: Int,
        treatmentCoverage: Coverage,
        treatmentScale: Double,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double
    ): Double {
        val chromosomeRange = ChromosomeRange(start, end, chromosome)
        val tc = treatmentCoverage.getBothStrandsCoverage(chromosomeRange) * treatmentScale
        return if (controlCoverage != null && controlScale != null) {
            val cc = controlCoverage.getBothStrandsCoverage(chromosomeRange) * controlScale
            max(0.0, tc - beta * cc)
        } else {
            tc
        }
    }

}

fun DoubleArray.standardDeviation(): Double {
    var sum = 0.0
    var sumSq = 0.0
    for (value in this) {
        sum += value
        sumSq += value * value
    }
    return sqrt((sumSq - sum * sum / size) / size)
}
