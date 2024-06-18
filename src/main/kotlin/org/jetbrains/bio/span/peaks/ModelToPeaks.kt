package org.jetbrains.bio.span.peaks

import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_STEPS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BACKGROUND_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BACKGROUND_SENSITIVITY_FRAGMENTATION_MILD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BACKGROUND_SENSITIVITY_FRAGMENTATION_MODERATE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BACKGROUND_SENSITIVITY_FRAGMENTATION_SEVERE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP_FRAGMENTATION_MILD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP_FRAGMENTATION_MODERATE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP_FRAGMENTATION_SEVERE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_LENGTH_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SIGNAL_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_MILD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_MODERATE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_PIVOT_THRESHOLD_SEVERE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanFitResults.Companion.LOG
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_GAPS_VARIANTS
import org.jetbrains.bio.span.statistics.util.PoissonUtil
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import kotlin.math.ceil
import kotlin.math.ln
import kotlin.math.ln1p
import kotlin.math.min


private const val COVERAGE_PERCENTILE_MIN = 5.0

private const val COVERAGE_PERCENTILE_MAX = 95.0

object ModelToPeaks {

    /**
     * Main method to compute peaks from the model.
     *
     * - Estimate posterior probabilities
     * - Merge candidate background bins with relaxed posterior error probability into blocks
     *    This mitigates the problem of wide marks peaks split on strong fdrs.
     * - Pick those blocks with at least one confident foreground stricter bin inside
     * - Assign p-value to each peak based on combined p-values for cores (consequent foreground bins).
     *    In case when control track is present, we use Poisson CDF to estimate log P-value;
     *    otherwise, an average log PEP (posterior error probability) for bins in blocks is used.
     *    N% top significant blocks scores are aggregated using length-weighted average as P for peak.
     * - Compute qvalues by peaks p-values, filter by alpha.
     * - Optional clipping to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    fun getPeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        bgSensitivity: Double?,
        gap: Double?,
        cancellableState: CancellableState? = null,
    ): SpanPeaksResult {
        val fitInfo = spanFitResults.fitInfo
        // Prepare fit information for scores computations
        fitInfo.prepareData()

        val (bgSens2use, gap2use) =
            estimateParams(genomeQuery, spanFitResults, fdr, bgSensitivity, gap, cancellableState)

        // Collect candidates from model
        val candidates = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>() to emptyList()
            }
            getChromosomeCandidates(spanFitResults, chromosome, fdr, bgSens2use, gap2use)
        }

        // Estimate signal and noise average signal by candidates
        val canEstimateSignalToNoise = fitInfo is SpanAnalyzeFitInformation &&
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
        return SpanPeaksResult(fdr, bgSens2use, gap2use, peaks)
    }


    private fun estimateParams(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        fdr: Double,
        bgSensitivity: Double?,
        gap: Double?,
        cancellableState: CancellableState?
    ): Pair<Double, Double> {
        if (bgSensitivity != null && gap != null) {
            return bgSensitivity to gap
        }
        LOG.info("Estimating background sensitivity and gap...")
        val minPivotGap = SPAN_GAPS_VARIANTS.sorted().firstOrNull { g ->
            val (n1, al1, _) =
                estimateCandidatesNumberLenDist(genomeQuery, spanFitResults, fdr, 1.0, g, cancellableState)
            val (n01, al01, _) =
                estimateCandidatesNumberLenDist(genomeQuery, spanFitResults, fdr, 0.1, g, cancellableState)
            val (n0001, al0001, _) =
                estimateCandidatesNumberLenDist(genomeQuery, spanFitResults, fdr, 0.001, g, cancellableState)
            // Here we estimate the direction between sensitivity = 1.0, 0.1 and 0.001
            // Experimentally we observe that bad quality highly fragmented tracks has counterclockwise rotation
            // Counterclockwise direction can be checked as a positive sign of convex hull square formula
            val s = triangleSignedSquare(
                n1.toDouble(), al1,
                n01.toDouble(), al01,
                n0001.toDouble(), al0001
            )
            val sLn = triangleSignedSquare(
                ln1p(n1.toDouble()), ln1p(al1),
                ln1p(n01.toDouble()), ln1p(al01),
                ln1p(n0001.toDouble()), ln1p(al0001)
            )
            return@firstOrNull s > 0 || sLn > 0
        }
        LOG.info("Minimal pivot gap: $minPivotGap")
        var bgSensitivity2use = bgSensitivity ?: SPAN_DEFAULT_BACKGROUND_SENSITIVITY
        var gap2use = gap ?: SPAN_DEFAULT_GAP
        if (minPivotGap != null) {
            when {
                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_SEVERE -> {
                    LOG.warn("Severe peaks fragmentation detected.")
                    bgSensitivity2use = bgSensitivity ?: SPAN_DEFAULT_BACKGROUND_SENSITIVITY_FRAGMENTATION_SEVERE
                    gap2use = gap ?: SPAN_DEFAULT_GAP_FRAGMENTATION_SEVERE
                }

                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_MODERATE -> {
                    LOG.warn("Moderate peaks fragmentation detected.")
                    bgSensitivity2use = bgSensitivity ?: SPAN_DEFAULT_BACKGROUND_SENSITIVITY_FRAGMENTATION_MODERATE
                    gap2use = gap ?: SPAN_DEFAULT_GAP_FRAGMENTATION_MODERATE
                }

                minPivotGap <= SPAN_GAP_PIVOT_THRESHOLD_MILD -> {
                    LOG.warn("Mild peaks fragmentation detected.")
                    bgSensitivity2use = bgSensitivity ?: SPAN_DEFAULT_BACKGROUND_SENSITIVITY_FRAGMENTATION_MILD
                    gap2use = gap ?: SPAN_DEFAULT_GAP_FRAGMENTATION_MILD
                }

                else -> LOG.warn("Barely fragmentation detected. Using defaults.")
            }
        }
        LOG.info("Sensitivity: $bgSensitivity2use Gap: $gap2use")
        return bgSensitivity2use to gap2use
    }

    private fun triangleSignedSquare(x1: Double, y1: Double, x2: Double, y2: Double, x3: Double, y3: Double) =
        x1 * y2 - x2 * y1 + x2 * y3 - x3 * y2 + x3 * y1 - x1 * y3

    fun estimateCandidatesNumberLenDist(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        fdr: Double,
        backgroundSensitivity: Double,
        gap: Double,
        cancellableState: CancellableState?
    ): Triple<Int, Double, Double> {
        val candidates = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            if (!spanFitResults.fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>() to emptyList()
            }
            getChromosomeCandidates(spanFitResults, chromosome, fdr, backgroundSensitivity, gap)
        }
        val candidatesN = genomeQuery.get().sumOf { candidates[it].first.size }
        val candidatesL = genomeQuery.get().sumOf {
            candidates[it].first.sumOf { r ->
                r.length().toLong()
            }
        }
        val candidatesAL = if (candidatesN > 0) candidatesL.toDouble() / candidatesN else 0.0
        val distances = DoubleArray(candidatesN)
        var i = 0
        genomeQuery.get().forEach { chr ->
            var prevEnd = -1
            candidates[chr].first.forEach { (start, end) ->
                if (prevEnd != -1) {
                    distances[i++] = (start - prevEnd).toDouble()
                }
                prevEnd = end
            }
        }
        val pLow = StatUtils.percentile(distances, COVERAGE_PERCENTILE_MIN)
        val pHigh = StatUtils.percentile(distances, COVERAGE_PERCENTILE_MAX)
        var sumP = 0.0
        i = 0
        distances.forEach {
            if (it in pLow..pHigh) {
                sumP += it
                i += 1
            }
        }
        val candidatesAD = if (i != 0) sumP / i else 0.0
        return Triple(candidatesN, candidatesAL, candidatesAD)
    }


    fun getChromosomeCandidates(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome,
        fdr: Double,
        bgSensitivity: Double,
        gap: Double,
        scoreBlocks: Double = SPAN_SCORE_BLOCKS,
    ): Pair<List<Range>, List<List<Range>>> {
        require(bgSensitivity >= 0) { "Sensitivity should be >=0, got $bgSensitivity" }
        require(scoreBlocks in 0.0..1.0) { "Blocks for score computation should be in range 0-1, got $scoreBlocks" }

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
        val candidateBins = BitList(logNullMemberships.size) { logNullMemberships[it] <= logFdr * bgSensitivity }

        // Background candidates with foreground inside
        var candidates = candidateBins.aggregate()
        if (candidates.isEmpty()) {
            return emptyList<Range>() to emptyList()
        }

        // merge candidates by min relative distance
        var lastEnd = -1
        var lastD = Int.MAX_VALUE
        for ((start, end) in candidates) {
            val d = ceil((end - start) * gap).toInt()
            if (lastEnd != -1) {
                if (start - lastEnd < min(lastD, d)) {
                    candidateBins.set(lastEnd, start)
                }
            }
            lastEnd = end
            lastD = d
        }
        candidates = candidateBins.aggregate()

        val candidateBlocks = candidates.map { (start, end) ->
            val p = StatUtils.percentile(
                logNullMemberships.slice(start, end).toDoubleArray(), scoreBlocks * 100
            )
            val currentPeakBlocks = BitList(end - start) {
                logNullMemberships[it + start] <= p
            }.aggregate(0).map { Range(start + it.startOffset, start + it.endOffset) }
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
                val (cs, ce) = clipPeakByScore(chromosome, start, end, fitInfo, avgSignalDensity, avgNoiseDensity)
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