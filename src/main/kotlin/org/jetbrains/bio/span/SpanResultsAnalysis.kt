package org.jetbrains.bio.span

import com.google.common.collect.MinMaxPriorityQueue
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.FixedFragment
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.AUTOCORRELATION_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_BAD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_PROBLEMATIC
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.peaks.ModelToPeaks.adjustQualitySensitivityGap
import org.jetbrains.bio.span.peaks.ModelToPeaks.computeCorrelations
import org.jetbrains.bio.span.peaks.ModelToPeaks.detectGapModel
import org.jetbrains.bio.span.peaks.ModelToPeaks.detectSensitivity
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateCandidatesNumberLen
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateMinPivotGap
import org.jetbrains.bio.span.peaks.ModelToPeaks.getLogNullPvals
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_BACKGROUND_SENSITIVITY_VARIANTS
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_GAPS_VARIANTS
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.io.BufferedWriter
import java.io.FileWriter
import java.nio.file.Path
import kotlin.math.ceil
import kotlin.math.floor
import kotlin.math.max
import kotlin.math.sqrt

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

        LOG.info("Analysing log null pvalues distribution")
        val logNullPvals = getLogNullPvals(genomeQuery, spanFitResults, blacklistPath)
        logInfo("LogNullPVals std: ${logNullPvals.standardDeviation()}", infoWriter)

        prepareLogNullsTsvFile(logNullPvals, peaksPath)

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
        logInfo("Track roughness score: ${"%.3f".format(roughness)}", infoWriter)
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
        val sensitivityInfo = detectSensitivity(genomeQuery, spanFitResults, fdr)
        val (detailedSensitivities, candNs, candALs, t1, t2, t3, area) = sensitivityInfo
        LOG.info("t1 $t1, t2 $t2, t3 $t3")
        logInfo("Sensitivity point 1 ${detailedSensitivities[t1]}", infoWriter)
        logInfo("Sensitivity point 2 ${detailedSensitivities[t2]}", infoWriter)
        logInfo("Sensitivity point 3 ${detailedSensitivities[t3]}", infoWriter)
        logInfo("Sensitivity triangle area $area", infoWriter)

        val minPivotGap = estimateMinPivotGap(genomeQuery, spanFitResults, fdr, sensitivityInfo)
        logInfo("Minimal pivot gap: $minPivotGap", infoWriter)
        if (minPivotGap != null) {
            when {
                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_BAD ->
                    logInfo("Quality: bad", infoWriter)

                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_PROBLEMATIC ->
                    logInfo("Quality: problematic", infoWriter)

                else ->
                    logInfo("Quality: good", infoWriter)
            }
        } else {
            logInfo("Quality: good", infoWriter)
        }

        val estimatedSensitivity = sensitivityInfo.sensitivities[sensitivityInfo.t2]
        val (sensitivity2use, gap2use) = adjustQualitySensitivityGap(
            null,
            null,
            minPivotGap,
            estimatedSensitivity,
            maxGapModel
        )
        logInfo("Estimated sensitivity: $estimatedSensitivity", infoWriter)
        logInfo("Estimated gap: $maxGapModel", infoWriter)
        logInfo("Adjusted sensitivity: $sensitivity2use", infoWriter)
        logInfo("Adjusted gap: $gap2use", infoWriter)

        infoWriter?.close()

        prepareSensitivitiesTsvFile(genomeQuery, spanFitResults, peaksPath, detailedSensitivities, candNs, candALs, fdr)

        prepareSegmentsFiles(genomeQuery, spanFitResults, sensitivityInfo, peaksPath, fdr)
    }

    private fun prepareLogNullsTsvFile(
        logNullPvals: DoubleArray,
        peaksPath: Path?
    ) {
        LOG.info("Analysing log nulls percentiles...")
        logNullPvals.sort()
        val logNullPsFile = if (peaksPath != null) "$peaksPath.logps.tsv" else null
        logNullPsFile?.toPath()?.deleteIfExists()
        if (logNullPsFile != null) {
            LOG.info("See $logNullPsFile")
        }

        val logNullPsWriter = if (logNullPsFile != null)
            BufferedWriter(FileWriter(logNullPsFile))
        else
            null
        logInfo("Q\tLogNullP", logNullPsWriter, false)
        var q = 0.0
        var step = 1e-5
        while (1 / step > logNullPvals.size) {
            step *= 10
        }
        while (q < 0.005) {
            logInfo("$q\t${logNullPvals[(logNullPvals.size * q).toInt()]}", logNullPsWriter, false)
            q += step
        }
        logNullPsWriter?.close()
    }

    private fun prepareSensitivitiesTsvFile(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        peaksPath: Path?,
        detailedSensitivities: DoubleArray,
        candNs: IntArray,
        candALs: DoubleArray,
        fdr: Double
    ) {
        LOG.info("Analysing candidates characteristics wrt sensitivity and gap...")
        val sensDetailsFile = if (peaksPath != null) "$peaksPath.sensitivity.tsv" else null
        if (sensDetailsFile != null) {
            LOG.info("See $sensDetailsFile")
        }
        sensDetailsFile?.toPath()?.deleteIfExists()
        val sensDetailsWriter = if (sensDetailsFile != null)
            BufferedWriter(FileWriter(sensDetailsFile))
        else
            null
        logInfo("Sensitivity\tGap\tCandidatesN\tCandidatesAL", sensDetailsWriter, false)
        // Print extended sensitivity table for gap = 0
        for ((i, s) in detailedSensitivities.withIndex()) {
            val candidatesN = candNs[i]
            val candidatesAL = candALs[i]
            logInfo("$s\t0\t$candidatesN\t$candidatesAL", sensDetailsWriter, false)
        }
        for (g in SPAN_GAPS_VARIANTS) {
            if (g == 0) {
                continue
            }
            for (s in SPAN_BACKGROUND_SENSITIVITY_VARIANTS) {
                val (candidatesN, candidatesAL) =
                    estimateCandidatesNumberLen(genomeQuery, spanFitResults, fdr, s, g, null)
                logInfo("$s\t$g\t$candidatesN\t$candidatesAL", sensDetailsWriter, false)
            }
        }
        sensDetailsWriter?.close()
    }

    private fun prepareSegmentsFiles(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        sensitivityInfo: ModelToPeaks.SensitivityInfo,
        peaksPath: Path?,
        fdr: Double,
    ) {
        if (peaksPath == null) {
            return
        }
        LOG.info("Analysing sensitivity segments...")
        val t0 = 0
        val t1 = sensitivityInfo.t1
        val t2 = sensitivityInfo.t2
        val t3 = sensitivityInfo.t3
        val t4 = sensitivityInfo.sensitivities.size - 1
        val segments = intArrayOf(t4, t3, t2, t1, t0)
        segments.forEachIndexed { i, t ->
            val s = "%.2e".format(sensitivityInfo.sensitivities[t])
            val segmentPath = "$peaksPath.segment${i}_$s.bed"
            LOG.info("Preparing $segmentPath")
            val peaks = ModelToPeaks.getPeaks(
                spanFitResults, genomeQuery, fdr, sensitivityInfo.sensitivities[t], 0, false,
                CancellableState.current()
            ).toList()
            Peak.savePeaks(
                peaks, segmentPath.toPath(),
                "segment${i}_$s"
            )
        }
    }

    private fun logInfo(msg: String, infoWriter: BufferedWriter?, useLog: Boolean = true) {
        if (useLog)
            SpanCLA.LOG.info(msg)
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
