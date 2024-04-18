package org.jetbrains.bio.span.peaks

import org.apache.commons.math3.stat.StatUtils
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_STEPS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MAXIMUM_COVERAGE_DELTA
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS_DISTANCE
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanFitResults.Companion.LOG
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.statistics.util.PoissonUtil
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.ln
import kotlin.math.min


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
    fun computeChromosomePeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        bgSensitivity: Double,
        clip: Double,
        cancellableState: CancellableState? = null
    ): List<Peak> {
        val fitInfo = spanFitResults.fitInfo
        fitInfo.prepareData()
        // Collect peaks from candidates
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            computeChromosomePeaks(spanFitResults, chromosome, fdr, bgSensitivity, clip)
        }
        return genomeQuery.get().flatMap { map[it] }
    }


    fun computeChromosomePeaks(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome,
        fdr: Double,
        bgSensitivity: Double,
        clip: Double,
        cancellableState: CancellableState? = null
    ): List<Peak> {
        // Check that we have information for requested chromosome
        val fitInfo = spanFitResults.fitInfo
        if (chromosome.name !in fitInfo.chromosomesSizes) {
            LOG.warn("Ignore ${chromosome.name}: model doesn't contain information")
            return emptyList()
        }

        fitInfo.prepareData()

        // Check genome contigs either by coverage info or by name
        // When contig's coverage is much bigger that genome-wide, SPAN calls too long peaks
        if (fitInfo is SpanAnalyzeFitInformation &&
            fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        ) {
            // We assume that the order of chromosomes is the same as in chromosome.sizes file
            val regularChromosome = chromosome.genome.toQuery().get().first()
            val mean = fitInfo.score(chromosome.chromosomeRange) / chromosome.length
            val meanRegular = fitInfo.score(regularChromosome.chromosomeRange) / regularChromosome.length
            if (abs(mean - meanRegular) > SPAN_MAXIMUM_COVERAGE_DELTA) {
                LOG.warn("Ignore ${chromosome.name}: coverage is different $mean vs $regularChromosome :$meanRegular")
                return emptyList()
            }
        } else if ('_' in chromosome.name || "random" in chromosome.name) {
            LOG.warn("Ignore ${chromosome.name}: chromosome name looks like contig")
            return emptyList()
        }

        val logNullMemberships =
            spanFitResults.logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)
        val offsets = fitInfo.offsets(chromosome)

        val candidates = computeCandidatePeaks(
            logNullMemberships, ln(fdr),
            bgSensitivity = bgSensitivity,
            scoreBlocks = SPAN_SCORE_BLOCKS
        )

        return peaksFromCandidates(
            candidates,
            chromosome,
            fitInfo,
            logNullMemberships,
            offsets,
            fdr,
            clip,
            cancellableState
        )
    }

    private fun peaksFromCandidates(
        candidates: Pair<List<Range>, List<List<Range>>>?,
        chromosome: Chromosome,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        fdr: Double,
        clip: Double,
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

        // Additionally, clip peaks by local coverage signal
        val doClip = clip > 0 &&
                fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        val (avgSignal, avgNoise) = if (doClip)
            estimateSignalAndNoiseDensity(
                chromosome, fitInfo, candidatePeaks, offsets
            ) else
            0.0 to 0.0
        var clipStart = 0L
        var clipEnd = 0L
        val maxClippedThreshold = clip
        val maxClippedFraction = clip / 2
        val maxClippedScore = avgNoise + maxClippedThreshold * (avgSignal - avgNoise)
        if (clip > 0) {
            LOG.debug("Signal density $avgSignal, noise density $avgNoise")
        }

        // Filter result peaks by Q values
        val peaksLogQValues = Fdr.qvalidate(peaksLogPvalues, logResults = true)
        val lnFdr = ln(fdr)
        val resultPeaks = candidatePeaks.mapIndexedNotNull { i, (from, to) ->
            cancellableState?.checkCanceled()
            val logPValue = peaksLogPvalues[i]
            val logQValue = peaksLogQValues[i]
            if (logPValue > lnFdr || logQValue > lnFdr) {
                return@mapIndexedNotNull null
            }
            val start = offsets[from]
            val end = if (to < offsets.size) offsets[to] else chromosome.length
            val (clippedStart, clippedEnd) = if (doClip)
                clipPeakByScore(chromosome, start, end, fitInfo, maxClippedScore, maxClippedFraction)
            else
                start to end
            clipStart += clippedStart - start
            clipEnd += end - clippedEnd
            Peak(
                chromosome = chromosome,
                startOffset = clippedStart,
                endOffset = clippedEnd,
                mlogpvalue = -logPValue / LOG_10,
                mlogqvalue = -logQValue / LOG_10,
                // Value is either coverage of fold change
                value = fitInfo.score(ChromosomeRange(clippedStart, clippedEnd, chromosome)),
                // Score should be proportional original q-value
                score = min(1000.0, -logQValue / LOG_10).toInt()
            )
        }
        return resultPeaks
    }

    fun computeCandidatePeaks(
        logNullMemberships: F64Array,
        logFdr: Double,
        bgSensitivity: Double,
        scoreBlocks: Double,
        minDistanceBetweenBlocks: Double = SPAN_SCORE_BLOCKS_DISTANCE
    ): Pair<List<Range>, List<List<Range>>> {
        require(bgSensitivity >= 0) { "Sensitivity should be >=0, got $bgSensitivity" }
        require(scoreBlocks in 0.0..1.0) { "Blocks for score computation should be in range 0-1, got $scoreBlocks" }
        val candidateBins = BitList(logNullMemberships.size) { logNullMemberships[it] <= logFdr * bgSensitivity }

        // Background candidates with foreground inside
        val candidates = candidateBins.aggregate()
        if (candidates.isEmpty()) {
            return emptyList<Range>() to emptyList()
        }

        val candidatePeaks = arrayListOf<Range>()
        val candidatePeaksBlocks = arrayListOf<List<Range>>()

        var currentPeakStart = -1
        var currentPeakEnd = -1
        for ((start, end) in candidates) {
            when {
                currentPeakStart == -1 -> {
                    currentPeakStart = start
                    currentPeakEnd = end
                }

                start < currentPeakEnd -> {
                    currentPeakEnd = end
                }

                else -> {
                    val currentPeak = Range(currentPeakStart, currentPeakEnd)
                    val minD = ceil(currentPeak.length() * minDistanceBetweenBlocks).toInt()
                    val p = StatUtils.percentile(
                        logNullMemberships.slice(currentPeakStart, currentPeakEnd).toDoubleArray(), scoreBlocks * 100
                    )
                    val currentPeakBlocks = BitList(currentPeakEnd - currentPeakStart) {
                        logNullMemberships[it + currentPeakStart] <= p
                    }.aggregate(minD)
                        .map { Range(currentPeakStart + it.startOffset, currentPeakStart + it.endOffset) }
                    candidatePeaks.add(currentPeak)
                    candidatePeaksBlocks.add(currentPeakBlocks)
                    currentPeakStart = start
                    currentPeakEnd = end
                }
            }
        }
        // Process last candidate peak
        val currentPeak = Range(currentPeakStart, currentPeakEnd)
        val p = StatUtils.percentile(
            logNullMemberships.slice(currentPeakStart, currentPeakEnd).toDoubleArray(), scoreBlocks * 100
        )
        val minD = ceil(currentPeak.length() * minDistanceBetweenBlocks).toInt()
        val currentPeakBlocks = BitList(currentPeakEnd - currentPeakStart) {
            logNullMemberships[it + currentPeakStart] <= p
        }.aggregate(minD)
            .map { Range(currentPeakStart + it.startOffset, currentPeakStart + it.endOffset) }
        candidatePeaks.add(currentPeak)
        candidatePeaksBlocks.add(currentPeakBlocks)

        return candidatePeaks to candidatePeaksBlocks
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
            // No significant bin within a candidate
            if (blocks.isEmpty()) {
                blocks = listOf(candidatePeaks[idx])
            }
            val (from, to) = candidatePeaks[idx]
            val start = offsets[from]
            val end = if (to < offsets.size) offsets[to] else chromosome.length
            val fullRange = ChromosomeRange(start, end, chromosome)
            // This can be beneficial in case of low-coverage depth,
            // when single short blocks may produce insignificant enrichment
            val fullLogPs = if (isTreatmentAndControlAvailable) {
                // Estimate enrichment vs local coverage in control track
                val peakTreatment = fitInfo.scaledTreatmentCoverage(fullRange)
                val peakControl = fitInfo.scaledControlCoverage(fullRange)!!
                // Use +1 as a pseudo count to compute Poisson CDF
                PoissonUtil.logPoissonCdf(
                    ceil(peakTreatment).toInt() + 1, ceil(peakControl) + 1
                )
            } else {
                // Average posterior log error probability for block
                (from until to).sumOf { logNullMemberships[it] } / (to - from)
            }
            val blocksLogPs = blocks.map { (from, to) ->
                if (!isTreatmentAndControlAvailable) {
                    // Average posterior log error probability for block
                    return@map (from until to).sumOf<Int> { logNullMemberships[it] } / (to - from)
                }
                // Estimate enrichment vs local coverage in control track
                val blockStart = offsets[from]
                val blockEnd = if (to < offsets.size) offsets[to] else chromosome.length
                val blockRange = ChromosomeRange(blockStart, blockEnd, chromosome)
                val peakTreatment = fitInfo.scaledTreatmentCoverage(blockRange)
                val peakControl = fitInfo.scaledControlCoverage(blockRange)!!
                // Use +1 as a pseudo count to compute Poisson CDF
                return@map PoissonUtil.logPoissonCdf(
                    ceil(peakTreatment).toInt() + 1, ceil(peakControl) + 1
                )
            }
            val blocksAverageLogPs = lengthWeightedScores(blocks, blocksLogPs)
            // Take most significant between full peak and blocks estimation
            min(fullLogPs, blocksAverageLogPs)
        }
        return peaksLogPvalues
    }


    private fun estimateSignalAndNoiseDensity(
        chromosome: Chromosome,
        fitInfo: SpanFitInformation,
        peaks: List<Range>,
        offsets: IntArray
    ): Pair<Double, Double> {
        var sumSignalScore = 0.0
        var sumSignalLength = 0L
        var sumNoiseScore = 0.0
        var sumNoiseLength = 0L
        var prevNoiseStart = 0
        peaks.forEach { (from, to) ->
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
        val avgSignalDensity = if (sumSignalLength > 0) sumSignalScore / sumSignalLength else 0.0
        val avgNoiseDensity = if (sumNoiseLength > 0) sumNoiseScore / sumNoiseLength else 0.0
        if (sumSignalLength != 0L && sumNoiseLength != 0L && avgSignalDensity <= avgNoiseDensity) {
            LOG.warn("Average signal density $avgSignalDensity <= average noise density $avgNoiseDensity")
            return avgNoiseDensity to avgSignalDensity
        }
        return avgSignalDensity to avgNoiseDensity
    }

    private fun clipPeakByScore(
        chromosome: Chromosome,
        start: Int,
        end: Int,
        fitInfo: SpanFitInformation,
        maxClippedScore: Double,
        maxClippedFraction: Double
    ): Pair<Int, Int> {
        val maxClippedLength = maxClippedFraction * (end - start)
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