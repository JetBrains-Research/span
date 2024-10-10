package org.jetbrains.bio.span

import com.google.common.collect.MinMaxPriorityQueue
import com.google.common.util.concurrent.AtomicDouble
import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTATION_MAX_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MIN_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SENSITIVITY_N
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.ModelToPeaks.analyzeAdditiveCandidates
import org.jetbrains.bio.span.peaks.ModelToPeaks.computeCorrelations
import org.jetbrains.bio.span.peaks.ModelToPeaks.detectSensitivityTriangle
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateCandidatesNumberLens
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateGap
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateGenomeSignalNoiseAverage
import org.jetbrains.bio.span.peaks.ModelToPeaks.getChromosomeCandidates
import org.jetbrains.bio.span.peaks.ModelToPeaks.getLogNullPvals
import org.jetbrains.bio.span.peaks.ModelToPeaks.getLogNulls
import org.jetbrains.bio.span.peaks.ModelToPeaks.linSpace
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_GAPS_VARIANTS
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.util.await
import org.jetbrains.bio.util.deleteIfExists
import org.jetbrains.bio.util.toPath
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.io.BufferedWriter
import java.io.FileWriter
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.*

object SpanResultsAnalysis {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun doDeepAnalysis(
        actualModelPath: Path,
        spanFitResults: SpanFitResults,
        fitInfo: SpanAnalyzeFitInformation,
        genomeQuery: GenomeQuery,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        blackListPath: Path?,
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
        val modelLowSignalToNoise = spanFitResults.model is NB2ZHMM &&
                spanFitResults.model.outOfSignalToNoiseRatioRangeDown
        logInfo("Model low signal to noise: $modelLowSignalToNoise", infoWriter)

        LOG.info("Analysing auto correlations...")
        val ncq = fitInfo.normalizedCoverageQueries!!.first()
        val (controlScale, beta, minCorrelation) = ncq.coveragesNormalizedInfo
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
            logInfo("Control coverage: $controlTotal", infoWriter)
            logInfo("Control scale: $controlScale", infoWriter)
            logInfo("Beta: $beta", infoWriter)
            logInfo("Min control correlation: $minCorrelation", infoWriter)
        }
        val blackList = if (blackListPath != null) {
            LOG.info("Loading blacklist regions: $blackListPath")
            LocationsMergingList.load(genomeQuery, blackListPath)
        } else null
        LOG.info("Analysing coverage distribution...")
        var coverage = computeCoverageScores(
            genomeQuery,
            treatmentCoverage, controlCoverage, controlScale,
            beta, fitInfo.binSize, blackList
        )
        logInfo("Coverage >0 %: ${(100.0 * coverage.count { it > 0 } / coverage.size).toInt()}", infoWriter)
        coverage = coverage.filter { it > 0 }.toDoubleArray()
        logInfo("Coverage >0 max: ${coverage.maxOrNull()}", infoWriter)
        logInfo("Coverage >0 mean: ${coverage.average()}", infoWriter)
        logInfo("Coverage >0 median: ${StatUtils.percentile(coverage, 50.0)}", infoWriter)
        logInfo("Coverage >0 std: ${coverage.standardDeviation()}", infoWriter)

        LOG.info("Analysing log null pvalues distribution...")
        val logNullPvals = getLogNullPvals(genomeQuery, spanFitResults, blackList)
        logInfo("LogNullPVals mean: ${logNullPvals.average()}", infoWriter)
        logInfo("LogNullPVals std: ${logNullPvals.standardDeviation()}", infoWriter)

        LOG.info("Analysing tracks roughness...")
        val roughness = computeAverageRoughness(
            genomeQuery, treatmentCoverage, controlCoverage, controlScale, beta, blackList
        )
        logInfo("Track roughness: ${"%.3f".format(roughness)}", infoWriter)

        // Collect candidates from model
        val logNullMembershipsMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                // HACK, it shouldn't be used
                // Cannot be null because of GenomeMap and empty F64Array is not supported
                return@genomeMap F64Array.full(1, 0.0)
            }
            getLogNulls(spanFitResults, chromosome)
        }
        val bitList2reuseMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
            BitList(logNullMembershipsMap[chromosome].length)
        }


        LOG.debug("Analysing autocorrelation...")
        val coverageCorrelations = computeCorrelations(coverage)
        val logNullPValsCorrelations = computeCorrelations(logNullPvals)
        val avgAutoCorrelation = logNullPValsCorrelations.average()
        logInfo("Average autocorrelation score: $avgAutoCorrelation", infoWriter)

        LOG.info("Analysing sensitivity...")
        val minLogNull = genomeQuery.get().minOf { logNullMembershipsMap[it].min() }
        // Limit value due to floating point errors
        val maxLogNull = min(SPAN_MIN_SENSITIVITY, genomeQuery.get().maxOf { logNullMembershipsMap[it].max() })
        val sensitivities = linSpace(minLogNull, maxLogNull, SPAN_SENSITIVITY_N)
        val si = detectSensitivityTriangle(
            genomeQuery, spanFitResults, logNullMembershipsMap, bitList2reuseMap, sensitivities
        )
        val sensitivity2use: Double
        when {
            sensitivityCmdArg != null ->
                sensitivity2use = sensitivityCmdArg
            si != null -> {
                val (beforeMerge, stable, beforeNoise) = si
                logInfo("Sensitivity beforeMerge: ${sensitivities[beforeMerge]}", infoWriter)
                logInfo("Sensitivity beforeMerge index: $beforeMerge", infoWriter)
                logInfo("Sensitivity stable: ${sensitivities[stable]}", infoWriter)
                logInfo("Sensitivity stable index: $stable", infoWriter)
                logInfo("Sensitivity beforeNoise: ${sensitivities[beforeNoise]}", infoWriter)
                logInfo("Sensitivity beforeNoise index: $beforeNoise", infoWriter)

                val sensitivitiesLimited =
                    sensitivities.slice(si.beforeMerge until si.stable).toDoubleArray()
                val (totals, news) = analyzeAdditiveCandidates(
                    genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                    sensitivitiesLimited, true
                )
                val minAdditionalIdx = sensitivitiesLimited.indices
                    .minByOrNull { news[it].toDouble() / totals[it].toDouble() }!!
                val minAdditionalSensitivity = sensitivitiesLimited[minAdditionalIdx]
                logInfo("Minimal additional: $minAdditionalSensitivity", infoWriter)
                logInfo("Minimal additional index: ${si.beforeMerge + minAdditionalIdx}", infoWriter)
                sensitivity2use = minAdditionalSensitivity
            }
            else -> {
                LOG.error("Failed to automatically estimate sensitivity")
                sensitivity2use = SPAN_DEFAULT_SENSITIVITY
            }
        }
        logInfo("Sensitivity2use: $sensitivity2use", infoWriter)

        LOG.debug("Analysing fragmentation...")
        val candidateGapNs = IntArray(SPAN_FRAGMENTATION_MAX_GAP) {
            estimateCandidatesNumberLens(
                genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                sensitivity2use, it
            ).n
        }
        val avgFragmentationScore = candidateGapNs.sumOf { it.toDouble() / candidateGapNs[0] }
        logInfo("Average fragmentation score: $avgFragmentationScore", infoWriter)

        val gap2use = if (gapCmdArg != null) {
            gapCmdArg
        } else {
            estimateGap(avgAutoCorrelation, avgFragmentationScore)
        }

        logInfo("Gap2use: $gap2use", infoWriter)

        val candidatesMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(chromosome, logNullMemberships, bitList2reuse, sensitivity2use, gap2use)
        }

        val (avgSignalDensity, avgNoiseDensity) =
            estimateGenomeSignalNoiseAverage(genomeQuery, fitInfo, candidatesMap, true)
        logInfo("Candidates signal density: $avgSignalDensity", infoWriter)
        logInfo("Candidates noise density: $avgNoiseDensity", infoWriter)
        val signalToNoise = if (avgNoiseDensity != 0.0) avgSignalDensity / avgNoiseDensity else 0.0
        logInfo("Coverage signal to noise: $signalToNoise", infoWriter)

        if (fitInfo.isControlAvailable()) {
            val signalToControl = computeSignalToControlAverage(genomeQuery, fitInfo, candidatesMap, true)
            logInfo("Coverage signal to control: $signalToControl", infoWriter)
        }
        infoWriter?.close()

        prepareCoveragePercentilesTsvFile(coverage, peaksPath)

        prepareLogNullPercentilesTsvFile(logNullPvals, peaksPath)

        prepareAutocorrelationTsvFile(coverageCorrelations, ".ac.coverage.tsv", peaksPath)

        prepareAutocorrelationTsvFile(logNullPValsCorrelations, ".ac.pvals.tsv", peaksPath)

        prepareFragmentationTsvFile(peaksPath, candidateGapNs)

        prepareSensitivitiesTsvFile(
            genomeQuery, spanFitResults, logNullMembershipsMap, bitList2reuseMap,
            peaksPath, sensitivities
        )

        prepareSegmentsTsvFile(genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap, sensitivities, peaksPath)
    }

    private fun prepareAutocorrelationTsvFile(correlations: DoubleArray, suffix: String, peaksPath: Path?) {
        val autocorrelationFile = if (peaksPath != null) "$peaksPath$suffix" else null
        autocorrelationFile?.toPath()?.deleteIfExists()
        if (autocorrelationFile != null) {
            LOG.info("See $autocorrelationFile")
        }

        val autocorrelationWriter = if (autocorrelationFile != null)
            BufferedWriter(FileWriter(autocorrelationFile))
        else
            null
        logInfo("D\tCorrelation", autocorrelationWriter, false)
        correlations.forEachIndexed { i, d ->
            logInfo("${i + 1}\t$d", autocorrelationWriter, false)
        }
        autocorrelationWriter?.close()
    }

    private fun prepareCoveragePercentilesTsvFile(
        coverage: DoubleArray,
        peaksPath: Path?
    ) {
        LOG.info("Analysing coverage percentiles...")
        val coveragePercFile = if (peaksPath != null) "$peaksPath.coverage.tsv" else null
        coveragePercFile?.toPath()?.deleteIfExists()
        if (coveragePercFile != null) {
            LOG.info("See $coveragePercFile")
        }

        val coveragePercWriter = if (coveragePercFile != null)
            BufferedWriter(FileWriter(coveragePercFile))
        else
            null
        logInfo("Q\tCoverage", coveragePercWriter, false)
        for (i in 0 until 100) {
            logInfo("${i + 1}\t${StatUtils.percentile(coverage, i + 1.0)}", coveragePercWriter, false)
        }
        coveragePercWriter?.close()
    }

    private fun prepareLogNullPercentilesTsvFile(
        logNullPvals: DoubleArray,
        peaksPath: Path?,
        maxQ: Double = 0.01,
        step: Double = 1e-5,
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
        var realStep = step
        while (1 / step > logNullPvals.size) {
            realStep *= 10
        }
        var q = 0.0
        while (q < maxQ) {
            logInfo("$q\t${logNullPvals[(logNullPvals.size * q).toInt()]}", logNullPsWriter, false)
            q += realStep
        }
        logNullPsWriter?.close()
    }

    private fun prepareFragmentationTsvFile(
        peaksPath: Path?,
        candidateGapNs: IntArray,
    ) {
        LOG.info("Analysing candidates wrt gap...")
        val gapsDetailsFile = if (peaksPath != null) "$peaksPath.gaps.tsv" else null
        if (gapsDetailsFile != null) {
            LOG.info("See $gapsDetailsFile")
        }
        gapsDetailsFile?.toPath()?.deleteIfExists()
        val gapsDetailsWriter = if (gapsDetailsFile != null)
            BufferedWriter(FileWriter(gapsDetailsFile))
        else
            null
        logInfo(
            "Gap\tCandidatesN",
            gapsDetailsWriter, false
        )
        candidateGapNs.forEachIndexed { g, n ->
            logInfo("$g\t$n", gapsDetailsWriter, false)
        }
        gapsDetailsWriter?.close()
    }

    private fun prepareSensitivitiesTsvFile(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        peaksPath: Path?,
        sensitivities: DoubleArray,
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
        logInfo(
            "Sensitivity\tGap\tCandidatesN\tCandidatesAL\tCandidatesML\tSignalNoiseRatio\tSignalControlRatio",
            sensDetailsWriter, false
        )
        for (g in SPAN_GAPS_VARIANTS) {
            for (s in sensitivities) {
                val candidatesMap = genomeMap(genomeQuery, parallel = true) { chromosome ->
                    if (!spanFitResults.fitInfo.containsChromosomeInfo(chromosome)) {
                        return@genomeMap emptyList<Range>()
                    }
                    val logNullMemberships = logNullMembershipsMap[chromosome]
                    val bitList2reuse = bitList2reuseMap[chromosome]
                    getChromosomeCandidates(chromosome, logNullMemberships, bitList2reuse, s, g)
                }
                val candidatesList = genomeQuery.get().flatMap { chromosome ->
                    candidatesMap[chromosome].map {
                        Location(it.startOffset, it.endOffset, chromosome)
                    }
                }
                val total = candidatesList.size
                val lengths = DoubleArray(total) { candidatesList[it].length().toDouble() }
                val avgL = lengths.average()
                val medianL = StatUtils.percentile(lengths, 50.0)
                val signalToNoise: Double
                val signalToControl: Double
                if (g == 0) {
                    val (avgSignalDensity, avgNoiseDensity) =
                        estimateGenomeSignalNoiseAverage(genomeQuery, spanFitResults.fitInfo, candidatesMap, true)
                    signalToNoise = if (avgNoiseDensity != 0.0) avgSignalDensity / avgNoiseDensity else 0.0
                    signalToControl = if (spanFitResults.fitInfo.isControlAvailable()) {
                        computeSignalToControlAverage(genomeQuery, spanFitResults.fitInfo, candidatesMap, true)
                    } else 0.0
                } else {
                    signalToNoise = 0.0
                    signalToControl = 0.0
                }
                logInfo("$s\t$g\t$total\t$avgL\t$medianL\t$signalToNoise\t$signalToControl", sensDetailsWriter, false)
            }
        }
        sensDetailsWriter?.close()
    }

    private fun prepareSegmentsTsvFile(
        genomeQuery: GenomeQuery,
        spanFitInformation: SpanFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivities: DoubleArray,
        peaksPath: Path?,
    ) {
        val (totals, news) =
            analyzeAdditiveCandidates(
                genomeQuery, spanFitInformation, logNullMembershipsMap, bitList2reuseMap,
                sensitivities, true
            )

        LOG.info("Analysing segments...")
        val segmentsFile = if (peaksPath != null) "$peaksPath.segments.tsv" else null
        if (segmentsFile != null) {
            LOG.info("See $segmentsFile")
        }
        segmentsFile?.toPath()?.deleteIfExists()
        val segmentsWriter = if (segmentsFile != null)
            BufferedWriter(FileWriter(segmentsFile))
        else
            null
        logInfo(
            "Sensitivity\tCandidatesN\tNew\tOld",
            segmentsWriter, false
        )
        sensitivities.forEachIndexed { i, s ->
            val total = totals[i]
            val new = news[i]
            val old = total - new
            logInfo(
                "$s\t$total\t$new\t$old",
                segmentsWriter, false
            )
        }
        segmentsWriter?.close()
    }

    private fun logInfo(msg: String, infoWriter: BufferedWriter?, useLog: Boolean = true) {
        if (useLog)
            SpanCLA.LOG.info(msg)
        if (infoWriter != null) {
            infoWriter.write(msg)
            infoWriter.newLine()
        }
    }

    val REGION_LEN = 10_000
    val TOP_REGIONS = 10_000
    val WORK_REGIONS = 200
    val RESOLUTION = 100

    /**
     * Detects roughness of the coverage.
     */

    private fun computeAverageRoughness(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        blackList: LocationsMergingList? = null,
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
                    treatmentCoverage, controlCoverage, controlScale, beta
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
                    treatmentCoverage, controlCoverage, controlScale, beta
                )
            }
            val mean = stats.average()
            val std = stats.standardDeviation()
            std_means[i / step] = if (mean > 0) std / mean else 0.0
        }
        return std_means.average()
    }


    private fun computeSignalToControlAverage(
        genomeQuery: GenomeQuery,
        fitInfo: SpanFitInformation,
        candidates: GenomeMap<List<Range>>,
        parallel: Boolean
    ): Double {
        val canEstimateSignalToNoise = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        if (!canEstimateSignalToNoise) {
            return 0.0
        }
        val signalToControls = AtomicDouble()
        val signalToControlsN = AtomicInteger()
        // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
        // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
        val treatmentCovs = (fitInfo as SpanAnalyzeFitInformation)
            .normalizedCoverageQueries!!.map { it.treatmentReads.get() }
        val controlCovs = fitInfo.normalizedCoverageQueries!!.map { it.controlReads!!.get() }
        val controlScales = fitInfo.normalizedCoverageQueries!!.map { it.coveragesNormalizedInfo.controlScale }

        genomeQuery.get().map { chromosome ->
            Callable {
                if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                    return@Callable
                }
                val offsets = fitInfo.offsets(chromosome)
                candidates[chromosome].forEach { (from, to) ->
                    val start = offsets[from]
                    val end = if (to < offsets.size) offsets[to] else chromosome.length
                    val chromosomeRange = ChromosomeRange(start, end, chromosome)
                    // val score = fitInfo.score(chromosomeRange)
                    val score = SpanAnalyzeFitInformation.fastScore(treatmentCovs, chromosomeRange)
                    // val controlScore = fitInfo.controlScore(chromosomeRange)
                    val controlScore = SpanAnalyzeFitInformation.fastControlScore(controlCovs, controlScales, chromosomeRange)
                    val ratio = if (controlScore != 0.0) score / controlScore else 0.0
                    signalToControls.addAndGet(ratio)
                    signalToControlsN.addAndGet(1)
                }
            }
        }.await(parallel)
        return if (signalToControls.get() > 0) signalToControls.get() / signalToControls.get() else 0.0
    }


    private fun computeCoverageScores(
        genomeQuery: GenomeQuery,
        treatmentCoverage: Coverage,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double,
        bin: Int,
        blackList: LocationsMergingList?
    ): DoubleArray {
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
                    treatmentCoverage, controlCoverage, controlScale, beta
                )
            }
        }
        if (blackList != null) {
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, totalBins)
        }
        return coverage
    }


    private fun coverage(
        chromosome: Chromosome,
        start: Int, end: Int,
        treatmentCoverage: Coverage,
        controlCoverage: Coverage?,
        controlScale: Double?,
        beta: Double
    ): Double {
        val chromosomeRange = ChromosomeRange(start, end, chromosome)
        val tc = treatmentCoverage.getBothStrandsCoverage(chromosomeRange).toDouble()
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
