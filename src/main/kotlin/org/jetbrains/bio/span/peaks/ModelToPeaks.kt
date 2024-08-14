package org.jetbrains.bio.span.peaks

import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.span.SpanResultsAnalysis
import org.jetbrains.bio.span.SpanResultsAnalysis.getCandidates
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_CHECKPOINT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_MAX_SHIFT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_NARROW_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_BROAD_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_STEPS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_LENGTH_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SIGNAL_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTATION_CHECKPOINT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_FRAGMENTATION_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_TRIANGLE_N
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_TRIANGLE_RELAXED_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_TRIANGLE_STRICT_SENSITIVITY
import org.jetbrains.bio.span.statistics.hmm.NB2ZHMM
import org.jetbrains.bio.span.statistics.util.PoissonUtil
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.*
import kotlin.math.*


object ModelToPeaks {

    val LOG: Logger = LoggerFactory.getLogger(ModelToPeaks.javaClass)

    /**
     * Main method to compute peaks from the model.
     *
     * 1) Estimate posterior probabilities.
     * 2) Merge candidate bins with relaxed posterior error probability (sensitivity) into candidates.
     *  This mitigates the problem of wide marks peaks split on strong fdrs.
     * 3) Find stable point in candidates number vs average length plot.
     * 4) Compute log null probabilities' autocorrelation.
     * 5) Increase gap between adjacent enriched bins with stable sensitivity to analyze fragmentation.
     * 6) Adjust sensitivity and gap according to autocorrelation and fragmentation of data.
     * 7) Assign p-value to each peak based on combined p-values for cores (consequent foreground bins).
     *    In case when control track is present, we use Poisson CDF to estimate log P-value;
     *    otherwise, an average log PEP (posterior error probability) for bins in blocks is used.
     *    N% top significant blocks scores are aggregated using length-weighted average as P for peak.
     * 8) Compute qvalues by peaks p-values, filter by alpha.
     * 10) Optional clipping to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    fun getPeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        clip: Boolean = true,
        blackListPath: Path? = null,
        cancellableState: CancellableState? = null,
    ): SpanPeaksResult {
        val fitInfo = spanFitResults.fitInfo
        // Prepare fit information for scores computations
        fitInfo.prepareData()

        val (sensitivity2use, gap2use) = estimateParams(
            genomeQuery, spanFitResults, fdr, sensitivityCmdArg, gapCmdArg,
            blackListPath, cancellableState
        )

        // Collect candidates from model
        val candidates = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>() to emptyList()
            }
            getChromosomeCandidates(spanFitResults, chromosome, fdr, sensitivity2use, gap2use)
        }

        // Estimate signal and noise average signal by candidates
        val canEstimateSignalToNoise = clip && fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }

        val (avgSignalDensity, avgNoiseDensity) = if (canEstimateSignalToNoise)
            estimateGenomeSignalNoiseAverage(genomeQuery, fitInfo, candidates).apply {
                LOG.debug("Signal density $first, noise density $second")
            }
        else
            null to null


        // Collect peaks
        val peaks = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                return@genomeMap emptyList<Peak>()
            }
            val chrCandidates = candidates[chromosome]
            val logNullMemberships = getLogNulls(spanFitResults, chromosome)
            getChromosomePeaksFromCandidates(
                chromosome, chrCandidates, fitInfo,
                logNullMemberships, fitInfo.offsets(chromosome), fdr,
                avgSignalDensity, avgNoiseDensity,
                cancellableState = cancellableState
            )
        }
        return SpanPeaksResult(fdr, sensitivity2use, gap2use, peaks)
    }


    private fun estimateParams(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        fdr: Double,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        blackListPath: Path?,
        cancellableState: CancellableState?
    ): Pair<Double, Int> {
        if (sensitivityCmdArg != null && gapCmdArg != null) {
            return sensitivityCmdArg to gapCmdArg
        }
        cancellableState?.checkCanceled()

        val modelLowSignalToNoise = spanFitResults.model is NB2ZHMM &&
                spanFitResults.model.outOfSignalToNoiseRatioRangeDown
        LOG.info("Low signal to noise: $modelLowSignalToNoise")
        cancellableState?.checkCanceled()

        LOG.info("Analysing sensitivity and gap...")
        val sensitivityInfo = detectSensitivityTriangle(genomeQuery, spanFitResults, fdr)
        cancellableState?.checkCanceled()

        LOG.info("Analysing log null pvalues distribution...")
        val logNullPvals = getLogNullPvals(genomeQuery, spanFitResults, blackListPath)
        val autocorrelationScore = getAutocorrelation(logNullPvals, SPAN_AUTOCORRELATION_CHECKPOINT)
        LOG.info("Autocorrelation score: $autocorrelationScore")

        LOG.info("Analysing fragmentation...")
        val sensitivities = sensitivityInfo.sensitivities
        val s2 = sensitivities[sensitivityInfo.t2]
        val candidatesNGapI = getCandidates(genomeQuery, spanFitResults, fdr, s2, SPAN_FRAGMENTATION_CHECKPOINT).size
        val candidatesNGap0 = getCandidates(genomeQuery, spanFitResults, fdr, s2, 0).size
        val fragmentationScore = if (candidatesNGap0 == 0) 1.0 else candidatesNGapI.toDouble() / candidatesNGap0
        LOG.info("Fragmentation score: $fragmentationScore")

        val (sensitivity2use, gap2use) = actualSensitivityGap(
            sensitivityInfo, sensitivityCmdArg, gapCmdArg,
            autocorrelationScore, fragmentationScore
        )
        LOG.info("Actual sensitivity: $sensitivity2use gap: $gap2use")
        cancellableState?.checkCanceled()
        return sensitivity2use to gap2use
    }

    private fun triangleSignedSquare(
        x1: Double, y1: Double,
        x2: Double, y2: Double,
        x3: Double, y3: Double
    ) =
        x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3

    @Suppress("ArrayInDataClass")
    data class SensitivityInfo(
        val sensitivities: DoubleArray,
        val t1: Int, val t2: Int, val t3: Int, val triangleArea: Double,
    )

    fun detectSensitivityTriangle(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        fdr: Double,
        relaxed: Double = SPAN_TRIANGLE_RELAXED_SENSITIVITY,
        strict: Double = SPAN_TRIANGLE_STRICT_SENSITIVITY,
        n: Int = SPAN_TRIANGLE_N
    ): SensitivityInfo {
        LOG.debug("Compute sensitivity triangle...")
        val strictLog10 = log10(strict)
        val relaxedLog10 = log10(relaxed)
        val sensitivities =
            DoubleArray(n) { exp((relaxedLog10 + (strictLog10 - relaxedLog10) * it.toDouble() / (n - 1)) * ln(10.0)) }
        val candidatesNs = IntArray(n)
        val candidatesALs = DoubleArray(n)
        for ((i, s) in sensitivities.withIndex()) {
            val (candidatesN, candidatesAL) = estimateCandidatesNumberLen(genomeQuery, spanFitResults, fdr, s, 0)
            candidatesNs[i] = candidatesN
            candidatesALs[i] = candidatesAL
        }
        val im1 = distanceArgmin(sensitivities, 1e-4)
        val im2 = distanceArgmin(sensitivities, 10.0)
        LOG.debug("Limit im1: {}, im2: {}", im1, im2)
        var maxArea = 0.0
        var i3 = -1
        var i2 = -1
        var i1 = -1
        for (i in im1 until im2) {
            val i3spb = findSensitivityTriangleMaxAreaBetween(
                candidatesNs, candidatesALs, 0, i, -1
            )
            val i1spb = findSensitivityTriangleMaxAreaBetween(
                candidatesNs, candidatesALs, i, n - 1, -1
            )
            if (i3spb[1] == -1.0 || i1spb[1] == -1.0) {
                continue
            }
            // We want both parts to be balanced so geometric mean optimization is better here
            val area = sqrt(i1spb[3] * i3spb[3])
            if (area > maxArea) {
                maxArea = area
                i1 = i1spb[1].toInt()
                i3 = i3spb[1].toInt()
                i2 = i
            }
        }
        LOG.debug("Intermediate I3: {}, I2: {}, I1: {}", i3, i2, i1)
        if (i3 == -1 || i2 == -1 || i1 == -1) {
            LOG.warn("Failed to estimate sensitivity triangle")
            return SensitivityInfo(sensitivities, im2, (im1 + im2) / 2, im1, maxArea)
        }
        // Update i3, i1 points to be closer to i2 for more accurate pivot estimations
        val i3spb = findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, i3, i2, -1)
        if (i3spb[1].toInt() != -1) {
            i3 = i3spb[1].toInt()
        }
        val i1spb = findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, i2, i1, -1)
        if (i1spb[1].toInt() != -1) {
            i1 = i1spb[1].toInt()
        }
        val result = SensitivityInfo(sensitivities, i3, i2, i1, maxArea)
        LOG.debug(
            "Result I3: {}: {}, I2: {}: {}, I1: {}: {}",
            i3, sensitivities[i3], i2, sensitivities[i2], i1, sensitivities[i1]
        )
        return result
    }

    private fun findSensitivityTriangleMaxAreaBetween(
        candidatesNs: IntArray,
        candidatesALs: DoubleArray,
        i3: Int, i1: Int, sign: Int = 1
    ): DoubleArray {
        if (i3 > i1) {
            return findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, i1, i3, sign)
        }
        var maxI1 = -1
        var maxI2 = -1
        var maxI3 = -1
        var maxArea = 0.0

        val n1 = candidatesNs[i1]
        val al1 = candidatesALs[i1]
        val n3 = candidatesNs[i3]
        val al3 = candidatesALs[i3]

        for (i2 in i3 until i1) {
            val n2 = candidatesNs[i2]
            val al2 = candidatesALs[i2]
            var area2 = triangleSignedSquare(
                ln1p(n1.toDouble()), ln1p(al1),
                ln1p(n2.toDouble()), ln1p(al2),
                ln1p(n3.toDouble()), ln1p(al3)
            )
            if (area2 * sign > 0)
                continue
            area2 = abs(area2)
            if (area2 > maxArea) {
                maxI1 = i1
                maxI2 = i2
                maxI3 = i3
                maxArea = area2
            }
        }
        return doubleArrayOf(maxI1.toDouble(), maxI2.toDouble(), maxI3.toDouble(), maxArea)
    }

    private fun distanceArgmin(array: DoubleArray, x: Double, start: Int = -1, end: Int = -1): Int {
        var minD = Double.MAX_VALUE
        var minI = -1


        for (i in array.indices) {
            if (start != -1 && i < start || end != -1 && i >= end) {
                continue
            }
            val a = array[i]
            val d = abs(a - x)
            if (d < minD) {
                minD = d
                minI = i
            }
        }
        return minI
    }

    fun actualSensitivityGap(
        sensitivityInfos: SensitivityInfo,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        autocorrelationScore: Double,
        fragmentationScore: Double
    ): Pair<Double, Int> {
        val sensitivity2use = when {
            sensitivityCmdArg != null -> sensitivityCmdArg
            // Keep narrow signal narrow
            autocorrelationScore < SPAN_AUTOCORRELATION_NARROW_THRESHOLD ->
                sensitivityInfos.sensitivities[sensitivityInfos.t3]
            // Stable point of HMM
            else -> sensitivityInfos.sensitivities[sensitivityInfos.t2]
        }
        val gap2use = when {
            gapCmdArg != null -> gapCmdArg
            // Keep narrow signal narrow
            autocorrelationScore < SPAN_AUTOCORRELATION_NARROW_THRESHOLD -> 0
            // Minimal gap is required for fragmented broad mark
            else -> SPAN_DEFAULT_GAP +
                    max(0.0, (SPAN_GAP_FRAGMENTATION_THRESHOLD - fragmentationScore) * SPAN_BROAD_GAP).toInt()
        }
        return sensitivity2use to gap2use
    }

    fun analyzeAdditiveCandidates(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        sensitivityInfo: SensitivityInfo,
        fdr: Double,
    ): Pair<IntArray, IntArray> {
        // Collect candidates from model
        val sensitivities = sensitivityInfo.sensitivities
        val candidatesPrev = genomeQuery.get().flatMap { chromosome ->
            if (!spanFitResults.fitInfo.containsChromosomeInfo(chromosome)) {
                return@flatMap emptyList()
            }
            getChromosomeCandidates(
                spanFitResults, chromosome, fdr,
                sensitivities[sensitivities.size - 1], 0
            )
                .first.map {
                    Location(it.startOffset, it.endOffset, chromosome)
                }
        }

        val totals = IntArray(sensitivities.size)
        val news = IntArray(sensitivities.size)
        var llPrev = LocationsMergingList.create(genomeQuery, candidatesPrev)
        totals[sensitivities.size - 1] = candidatesPrev.size
        news[sensitivities.size - 1] = candidatesPrev.size
        var i = sensitivities.size - 2
        while (i >= 0) {
            val s = sensitivities[i]
            val candidates = genomeMap(genomeQuery, parallel = true) { chromosome ->
                if (!spanFitResults.fitInfo.containsChromosomeInfo(chromosome)) {
                    return@genomeMap emptyList<Range>() to emptyList()
                }
                getChromosomeCandidates(spanFitResults, chromosome, fdr, s, 0)
            }

            val candidatesList = genomeQuery.get().flatMap { chromosome ->
                if (!spanFitResults.fitInfo.containsChromosomeInfo(chromosome)) {
                    return@flatMap emptyList()
                }
                candidates[chromosome].first.map {
                    Location(it.startOffset, it.endOffset, chromosome)
                }
            }
            val total = candidatesList.size
            val old = candidatesList.count { llPrev.intersects(it) }
            val new = total - old
            totals[i] = total
            news[i] = new
            llPrev = LocationsMergingList.create(genomeQuery, candidatesList)
            i -= 1
        }
        return totals to news
    }

    fun getLogNullPvals(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        blackListPath: Path?
    ): DoubleArray {
        // Limit genome query to top non-empty chromosomes
        val chrs = genomeQuery.get()
            .filter { spanFitResults.logNullMemberships.containsKey(it.name) }
            .sortedByDescending { it.length }
            .take(3)
            .map { it.name }.toTypedArray()
        val limitedQuery = GenomeQuery(genomeQuery.genome, *chrs)
        val bin = spanFitResults.fitInfo.binSize
        val totalBins = limitedQuery.get()
            .sumOf { spanFitResults.logNullMemberships[it.name]!!.f64Array(SpanModelFitExperiment.NULL).size }
        val result = DoubleArray(totalBins) { 0.0 }
        var i = 0
        val blackList = if (blackListPath != null) LocationsMergingList.load(limitedQuery, blackListPath) else null
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            val logNullPeps = spanFitResults.logNullMemberships[chr.name]!!.f64Array(SpanModelFitExperiment.NULL)
            for (j in 0 until logNullPeps.size) {
                val start = j * bin
                val end = (j + 1) * bin
                // Ignore blackList regions
                if (blackList != null && blackList.intersects(Location(start, end, chr))) {
                    blackListIgnored++
                    continue
                }
                result[i++] = logNullPeps[j]
            }
        }
        if (blackList != null) {
            SpanResultsAnalysis.LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, totalBins)
        }
        return result
    }

    fun computeCorrelations(
        values: DoubleArray,
        maxCorrelationDistance: Int = SPAN_AUTOCORRELATION_MAX_SHIFT
    ): DoubleArray {
        val distanceCorrelations = DoubleArray(min(values.size / 2, maxCorrelationDistance) + 1)
        distanceCorrelations[0] = 1.0
        for (d in 1 until distanceCorrelations.size) {
            val correlation = getAutocorrelation(values, d)
            // Stop computing correlation after small one
            if (Precision.equals(correlation, 0.0, 1e-6)) {
                break
            }
            distanceCorrelations[d] = correlation
        }
        return distanceCorrelations
    }

    fun getAutocorrelation(values: DoubleArray, d: Int): Double {
        val original = DoubleArray(values.size - d - 1)
        System.arraycopy(values, 0, original, 0, values.size - d - 1)
        val shifted = DoubleArray(values.size - d - 1)
        System.arraycopy(values, d, shifted, 0, values.size - d - 1)
        Arrays.fill(shifted, values.size - d - 1, shifted.size, 0.0)
        var correlation = if (original.size >= 2)
            PearsonsCorrelation().correlation(original, shifted)
        else
            0.0
        // Safeguard against NaN from Pearson correlation for zero or small vectors,
        // See example at: https://github.com/JetBrains-Research/span/issues/34
        if (correlation.isNaN() || !correlation.isFinite()) {
            correlation = 0.0
        }
        return correlation
    }

    fun estimateCandidatesNumberLen(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        fdr: Double,
        sensitivity: Double,
        gap: Int
    ): Pair<Int, Double> {
        val candidates = getCandidates(genomeQuery, spanFitResults, fdr, sensitivity, gap)
        val candidatesN = candidates.size
        val candidatesL = candidates.sumOf { it.length().toLong() }
        val candidatesAL = if (candidatesN > 0) candidatesL.toDouble() / candidatesN else 0.0
        return candidatesN to candidatesAL
    }


    fun getChromosomeCandidates(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome,
        fdr: Double,
        sensitivity: Double,
        gap: Int,
    ): Pair<List<Range>, List<List<Range>>> {
        require(sensitivity >= 0) { "Sensitivity should be >=0, got $sensitivity" }
        require(gap >= 0) { "Gap should be >=0, got $gap" }

        // Check that we have information for requested chromosome
        val fitInfo = spanFitResults.fitInfo
        if (!fitInfo.containsChromosomeInfo(chromosome)) {
            LOG.trace("Ignore ${chromosome.name}: model doesn't contain information")
            return emptyList<Range>() to emptyList()
        }

        if ('_' in chromosome.name ||
            "random" in chromosome.name.lowercase() ||
            "un" in chromosome.name.lowercase()
        ) {
            LOG.trace("Ignore ${chromosome.name}: chromosome name looks like contig")
            return emptyList<Range>() to emptyList()
        }

        val logNullMemberships = getLogNulls(spanFitResults, chromosome)
        val logFdr = ln(fdr)
        val candidateBins = BitList(logNullMemberships.size) { logNullMemberships[it] <= logFdr * sensitivity }

        // Background candidates with foreground inside
        val candidates = candidateBins.aggregate(gap)
        if (candidates.isEmpty()) {
            return emptyList<Range>() to emptyList()
        }
        // Compute significant blocks within candidate
        val candidateBlocks = candidates.map { (start, end) ->
            val p = StatUtils.percentile(
                logNullMemberships.slice(start, end).toDoubleArray(), SPAN_SCORE_BLOCKS * 100
            )
            val currentPeakBlocks = BitList(end - start) {
                logNullMemberships[it + start] <= p
            }.aggregate(SPAN_DEFAULT_GAP).map { Range(start + it.startOffset, start + it.endOffset) }
            return@map currentPeakBlocks
        }
        return candidates to candidateBlocks
    }

    fun getLogNulls(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome
    ) = spanFitResults.logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)

    fun getChromosomePeaksFromCandidates(
        chromosome: Chromosome,
        candidates: Pair<List<Range>, List<List<Range>>>?,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        fdr: Double,
        avgSignalDensity: Double?,
        avgNoiseDensity: Double?,
        cancellableState: CancellableState?
    ): List<Peak> {
        // Compute candidate bins and peaks with relaxed background settings
        // 1) Return broad peaks in case of broad modifications even for strict FDR settings
        // 2) Mitigate the problem when the number of peaks for strict FDR is much bigger than for relaxed FDR
        val (candidatePeaks, candidatePeaksBlocks) = candidates!!
        if (candidatePeaks.isEmpty()) {
            return emptyList()
        }

        // We want two invariants from peaks pvalues:
        // 1) The stricter FDR, the fewer peaks with smaller average length
        // 2) Peaks should not disappear when relaxing FDR
        // Peak score is computed as length-weighted average p-value in its consequent enriched bins.
        val peaksLogPvalues = estimateCandidatesLogPs(
            chromosome,
            candidatePeaks,
            candidatePeaksBlocks,
            fitInfo,
            logNullMemberships,
            offsets,
            cancellableState
        )

        // Filter result peaks by Q values
        val peaksLogQValues = Fdr.qvalidate(peaksLogPvalues, logResults = true)
        val lnFdr = ln(fdr)
        val canEstimateScore = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        val resultPeaks = candidatePeaks.mapIndexedNotNull { i, (from, to) ->
            cancellableState?.checkCanceled()
            val logPValue = peaksLogPvalues[i]
            val logQValue = peaksLogQValues[i]
            if (logPValue > lnFdr || logQValue > lnFdr) {
                return@mapIndexedNotNull null
            }
            var start = offsets[from]
            var end = if (to < offsets.size) offsets[to] else chromosome.length
            if (canEstimateScore && avgSignalDensity != null && avgNoiseDensity != null) {
                val (cs, ce) = clipPeakByScore(
                    chromosome, start, end, fitInfo, avgSignalDensity, avgNoiseDensity
                )
                start = cs
                end = ce
            }
            // Value is either coverage of fold change
            val value = if (canEstimateScore)
                fitInfo.score(ChromosomeRange(start, end, chromosome))
            else
                -logQValue
            Peak(
                chromosome = chromosome,
                startOffset = start,
                endOffset = end,
                mlogpvalue = -logPValue / LOG_10,
                mlogqvalue = -logQValue / LOG_10,
                value = value,
                // Score should be proportional original q-value
                score = min(1000.0, -logQValue / LOG_10).toInt()
            )
        }
        return resultPeaks
    }

    private fun estimateCandidatesLogPs(
        chromosome: Chromosome,
        candidatePeaks: List<Range>,
        candidatePeaksBlocks: List<List<Range>>,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        cancellableState: CancellableState?
    ): F64Array {
        val isTreatmentAndControlAvailable = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.hasControlData() &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }

        val peaksLogPvalues = F64Array(candidatePeaks.size) { idx ->
            cancellableState?.checkCanceled()
            var blocks = candidatePeaksBlocks[idx]
            // Use the whole peaks for estimation may produce more significant results
            // beneficial for short peaks, TFs, ATAC-seq etc.
            if (blocks.size <= 1) {
                blocks = listOf(candidatePeaks[idx])
            }
            val blocksLogPs = blocks.map { (from, to) ->
                // Model posterior log error probability for block
                val modelLogPs = (from until to).sumOf { logNullMemberships[it] }
                if (!isTreatmentAndControlAvailable) {
                    return@map modelLogPs
                }
                val blockStart = offsets[from]
                val blockEnd = if (to < offsets.size) offsets[to] else chromosome.length
                val blockRange = ChromosomeRange(blockStart, blockEnd, chromosome)
                val peakTreatment = fitInfo.score(blockRange)
                val peakControl = fitInfo.controlScore(blockRange)
                // Use +1 as a pseudo count to compute Poisson CDF
                val signalLogPs = PoissonUtil.logPoissonCdf(
                    ceil(peakTreatment).toInt() + 1,
                    peakControl + 1
                )
                // Combine both model and signal estimations
                // These can be of different magnitude, geometric mean is more stable
//                return@map -sqrt(modelLogPs * signalLogPs)
                return@map (modelLogPs + signalLogPs) / 2
            }
            return@F64Array lengthWeightedScores(blocks, blocksLogPs)
        }
        return peaksLogPvalues
    }

    fun estimateGenomeSignalNoiseAverage(
        genomeQuery: GenomeQuery,
        fitInfo: SpanFitInformation,
        candidates: GenomeMap<Pair<List<Range>, List<List<Range>>>>,
    ): Pair<Double, Double> {
        val canEstimateSignalToNoise = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        if (!canEstimateSignalToNoise) {
            return 0.0 to 0.0
        }
        var sumSignalScore = 0.0
        var sumSignalLength = 0L
        var sumNoiseScore = 0.0
        var sumNoiseLength = 0L
        genomeQuery.get().forEach { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                return@forEach
            }
            val chrCandidates = candidates[chromosome]
            val offsets = fitInfo.offsets(chromosome)
            var prevNoiseStart = 0
            chrCandidates.first.forEach { (from, to) ->
                val start = offsets[from]
                val end = if (to < offsets.size) offsets[to] else chromosome.length
                sumSignalScore += fitInfo.score(ChromosomeRange(start, end, chromosome))
                sumSignalLength += end - start
                sumNoiseScore += fitInfo.score(ChromosomeRange(prevNoiseStart, start, chromosome))
                sumNoiseLength += start - prevNoiseStart
                prevNoiseStart = end
            }
            if (prevNoiseStart < chromosome.length) {
                sumNoiseScore += fitInfo.score(ChromosomeRange(prevNoiseStart, chromosome.length, chromosome))
                sumNoiseLength += chromosome.length - prevNoiseStart
            }
        }
        val avgSignalDensity = if (sumSignalLength > 0) sumSignalScore / sumSignalLength else 0.0
        val avgNoiseDensity = if (sumNoiseLength > 0) sumNoiseScore / sumNoiseLength else 0.0
        if (sumSignalLength != 0L && sumNoiseLength != 0L && avgSignalDensity <= avgNoiseDensity) {
            LOG.warn("Average signal density $avgSignalDensity <= average noise density $avgNoiseDensity")
        }
        return avgSignalDensity to avgNoiseDensity
    }


    private fun clipPeakByScore(
        chromosome: Chromosome,
        start: Int,
        end: Int,
        fitInfo: SpanFitInformation,
        avgSignalDensity: Double,
        avgNoiseDensity: Double,
        clipSignal: Double = SPAN_DEFAULT_SIGNAL_CLIP,
        clipLength: Double = SPAN_DEFAULT_LENGTH_CLIP,
    ): Pair<Int, Int> {
        if (avgSignalDensity <= avgNoiseDensity) {
            return start to end
        }
        // Additionally, clip peaks by local coverage signal
        val maxClippedScore = avgNoiseDensity + clipSignal * (avgSignalDensity - avgNoiseDensity)
        val maxClippedLength = (end - start) * (1 - clipLength) / 2

        // Try to change the left boundary
        val maxStart = start + maxClippedLength
        var clippedStart = start
        var step = SPAN_CLIP_STEPS.size - 1
        while (step >= 0 && clippedStart <= maxStart) {
            val newStart = clippedStart + SPAN_CLIP_STEPS[step]
            if (newStart > maxStart) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clipScore = fitInfo.score(ChromosomeRange(start, newStart, chromosome)) / (newStart - start)
            if (clipScore < maxClippedScore) {
                clippedStart = newStart
                step = min(step + 1, SPAN_CLIP_STEPS.size - 1)
            } else {
                step -= 1
            }
        }
        // Try to change the right boundary
        val minEnd = end - maxClippedLength
        var clippedEnd = end
        step = SPAN_CLIP_STEPS.size - 1
        while (step >= 0 && clippedEnd >= minEnd) {
            val newEnd = clippedEnd - SPAN_CLIP_STEPS[step]
            if (newEnd < minEnd) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clipScore = fitInfo.score(ChromosomeRange(newEnd, end, chromosome)) / (end - newEnd)
            if (clipScore < maxClippedScore) {
                clippedEnd = newEnd
                step = min(step + 1, SPAN_CLIP_STEPS.size - 1)
            } else {
                step -= 1
            }
        }
        return clippedStart to clippedEnd
    }


    /**
     * Summary length-weighted score for a peak.
     * Use length weighted mean to take into account the difference in block lengths
     */
    private fun lengthWeightedScores(
        blocks: List<Range>,
        scores: List<Double>,
    ): Double {
        require(blocks.size == scores.size) { "Different lengths of blocks and scores lists" }
        if (blocks.size == 1) {
            return scores.first()
        }
        val sum = KahanSum()
        var l = 0
        blocks.zip(scores).sortedBy { it.second }
            .forEach { (b, p) ->
                sum += p * (b.startOffset - b.endOffset)
                l += b.startOffset - b.endOffset
            }
        return sum.result() / l
    }

    private val LOG_10 = ln(10.0)

}