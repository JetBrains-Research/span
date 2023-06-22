package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Range
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanFitResults.Companion.LOG
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.statistics.util.PoissonUtil
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.ceil
import kotlin.math.ln
import kotlin.math.max
import kotlin.math.min


object ModelToPeaks {

    /**
     * Main method to compute peaks from model.
     *
     * 1) Estimate posterior probabilities
     * 2) Merge candidate bins with relaxed posterior error probability into blocks
     * This mitigates the problem of wide marks peaks split on strong fdrs.
     * 3) Pick those blocks with at least one confident enriched bin inside
     * 4) Using gap merge blocks into candidate peaks.
     * 5) Assign p-value to each peak based on combined p-values for cores(consequent strict enriched bins).
     *    In case when control track is present, we use Poisson CDF to estimate log P value,
     *    otherwise an average log PEP (posterior error probability) for bins in blocks is used.
     *    50% top significant blocks scores are aggregated using length-weighted average as P for peak.
     * 6) Compute qvalues by peaks p-values, filter by alpha.
     */
    fun computeChromosomePeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        gap: Int,
        clip: Boolean,
        cancellableState: CancellableState? = null
    ): List<Peak> {
        resetCounters()
        val progress = Progress { title = "Computing peaks fdr=$fdr gap=$gap" }
            .bounded(genomeQuery.get().size.toLong())
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            val chromosomePeaks =
                computeChromosomePeaks(
                    spanFitResults, chromosome, fdr, gap, clip,
                    cancellableState = cancellableState
                )
            progress.report(1)
            chromosomePeaks
        }
        progress.done()
        Peak.LOG.debug(
            "Total candidate bins/candidate/result peaks " +
                    "${candidateBinsCounter.get()}/${candidatePeaksCounter.get()}/${resultPeaksCounter.get()}"
        )
        return genomeQuery.get().flatMap { map[it] }
    }


    fun computeChromosomePeaks(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome,
        fdr: Double,
        gap: Int,
        clip: Boolean,
        cancellableState: CancellableState? = null
    ): List<Peak> {
        // Check that we have information for requested chromosome
        val chromosomePeaks = if (chromosome.name in spanFitResults.fitInfo.chromosomesSizes) {
            spanFitResults.fitInfo.prepareData()
            getChromosomePeaks(
                chromosome,
                spanFitResults.fitInfo,
                spanFitResults.logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
                spanFitResults.fitInfo.offsets(chromosome),
                fdr,
                gap,
                clip,
                cancellableState = cancellableState
            )
        } else {
            LOG.debug(
                "NO peaks information for chromosome: ${chromosome.name} " +
                        "in fitInfo ${spanFitResults.fitInfo.build}"
            )
            emptyList()
        }
        return chromosomePeaks
    }

    private fun getChromosomePeaks(
        chromosome: Chromosome,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        fdr: Double,
        gap: Int,
        clip: Boolean,
        cancellableState: CancellableState? = null
    ): List<Peak> {
        // Compute candidate bins and peaks with relaxed settings
        // Relaxed probability allows for:
        // 1) Return broad peaks in case of broad modifications even for strict FDR settings
        // 2) Mitigate the problem when number of peaks for strict FDR is much bigger than for relaxed FDR
        val logFdr = ln(fdr)
        val relaxedLogFdr = relaxedLogFdr(logFdr)
        require(relaxedLogFdr >= logFdr)
        val relaxedBins = BitList(logNullMemberships.size) { logNullMemberships[it] <= relaxedLogFdr }
        val strictBins = BitList(logNullMemberships.size) { logNullMemberships[it] <= logFdr }
        val (candidateBins, candidatePeaks, candidatePeaksCores) =
            computeBinsCoresAndPeaks(relaxedBins, strictBins, gap)
        if (candidatePeaks.isEmpty()) {
            return emptyList()
        }

        val isTreatmentAndControlAvailable = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.hasControlData() &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }

        // We want two invariants from peaks pvalues:
        // 1) The more strict FDR, the fewer peaks with smaller average length
        // 2) Peaks should not disappear when relaxing FDR
        // Peak score is computed as length-weighted average p-value in its consequent enriched bins.
        val peaksLogPvalues = F64Array(candidatePeaks.size) { idx ->
            cancellableState?.checkCanceled()
            var coreBlocks = candidatePeaksCores[idx]
            // No significant bin within candidate
            if (coreBlocks.isEmpty()) {
                coreBlocks = listOf(candidatePeaks[idx])
            }
            val blocksLogPs = coreBlocks.map { (from, to) ->
                cancellableState?.checkCanceled()
                val start = offsets[from]
                val end = if (to < offsets.size) offsets[to] else chromosome.length
                val chromosomeRange = ChromosomeRange(start, end, chromosome)
                if (isTreatmentAndControlAvailable) {
                    // Estimate enrichment vs local coverage in control track
                    val peakTreatment = fitInfo.scaledTreatmentCoverage(chromosomeRange)
                    val peakControl = fitInfo.scaledControlCoverage(chromosomeRange)!!
                    PoissonUtil.logPoissonCdf(
                        ceil(peakTreatment).toInt() + PSEUDO_COUNT, ceil(peakControl) + PSEUDO_COUNT
                    )
                } else {
                    // Fallback to average posterior log error probability for block
                    (from until to).sumOf { logNullMemberships[it] } / (to - from)
                }
            }
            lengthWeightedScores(coreBlocks, blocksLogPs)
        }

        // Additionally clip peaks by coverage
        val clipEnabled = clip &&
                fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        val (avgSignal, avgNoise) = if (clipEnabled)
            estimateSignalAndNoiseDensity(
                chromosome, fitInfo, candidatePeaks, offsets
            ) else
            0.0 to 0.0
        var totalClipStart = 0L
        var totalClipEnd = 0L
        val maxClippedScore = avgNoise + MAX_CLIPPED_DELTA * (avgSignal - avgNoise)
        if (clip) {
            LOG.debug("Signal density $avgSignal, noise density $avgNoise")
        }

        // Filter result peaks by Q values
        val peaksLogQValues = Fdr.qvalidate(peaksLogPvalues, logResults = true)
        val resultPeaks = candidatePeaks.mapIndexedNotNull { i, (from, to) ->
            cancellableState?.checkCanceled()
            val start = offsets[from]
            val end = if (to < offsets.size) offsets[to] else chromosome.length
            val logPValue = peaksLogPvalues[i]
            val logQValue = peaksLogQValues[i]
            if (logPValue > logFdr || logQValue > logFdr) {
                return@mapIndexedNotNull null
            }

            val (clippedStart, clippedEnd) = if (clipEnabled)
                clipPeakByScore(
                    chromosome, start, end, fitInfo,
                    (MAX_CLIPPED_LENGTH * (end - start)).toInt(), maxClippedScore
                )
            else start to end
            totalClipStart += clippedStart - start
            totalClipEnd += end - clippedEnd
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

        Peak.LOG.debug(
            "$chromosome: candidate bins / candidate/result peaks " +
                    "${candidateBins.cardinality()} / ${candidatePeaks.size} / ${resultPeaks.size};"

        )
        candidateBinsCounter.addAndGet(candidateBins.cardinality())
        candidatePeaksCounter.addAndGet(candidatePeaks.size)
        resultPeaksCounter.addAndGet(resultPeaks.size)
        return resultPeaks
    }

    fun computeBinsCoresAndPeaks(
        relaxedBins: BitList,
        strictBins: BitList,
        gap: Int
    ): Triple<BitList, List<Range>, List<List<Range>>> {
        require(relaxedBins.size() == strictBins.size()) { "Different size" }
        (0 until relaxedBins.size()).forEach {
            require(!(strictBins[it] && !relaxedBins[it])) { "Strict positions should be covered in relaxed" }
        }
        val cores = strictBins.aggregate()
        // Relaxed candidates with core inside
        val candidates = relaxedBins.aggregate().filter { (start, end) ->
            (start until end).any { strictBins[it] }
        }

        if (candidates.isEmpty()) {
            return Triple(relaxedBins, emptyList(), emptyList())
        }

        val candidatePeaks = arrayListOf<Range>()
        val candidatePeaksCores = arrayListOf<List<Range>>()
        var coreIdx = 0
        var currentPeakStart = -1
        var currentPeakEnd = -1
        var currentPeaksCores = arrayListOf<Range>()
        for ((start, end) in candidates) {
            when {
                currentPeakStart == -1 -> {
                    currentPeakStart = start
                    currentPeakEnd = end
                }

                start < currentPeakEnd + gap -> {
                    currentPeakEnd = end
                }

                else -> {
                    val currentPeak = Range(currentPeakStart, currentPeakEnd)
                    while (coreIdx < cores.size && currentPeak.contains(cores[coreIdx])) {
                        currentPeaksCores.add(cores[coreIdx])
                        coreIdx++
                    }
                    if (currentPeaksCores.isNotEmpty()) {
                        candidatePeaks.add(currentPeak)
                        candidatePeaksCores.add(currentPeaksCores)
                    }
                    currentPeakStart = start
                    currentPeakEnd = end
                    currentPeaksCores = arrayListOf()
                }
            }
        }
        // Process last candidate peak
        val currentPeak = Range(currentPeakStart, currentPeakEnd)
        while (coreIdx < cores.size && currentPeak.contains(cores[coreIdx])) {
            currentPeaksCores.add(cores[coreIdx])
            coreIdx++
        }
        if (currentPeaksCores.isNotEmpty()) {
            candidatePeaks.add(currentPeak)
            candidatePeaksCores.add(currentPeaksCores)
        }

        return Triple(relaxedBins, candidatePeaks, candidatePeaksCores)
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

    /**
     * Tries to reduce range by [CLIP_STEPS] from both sides while increasing score.
     */
    private fun clipPeakByScore(
        chromosome: Chromosome,
        start: Int,
        end: Int,
        fitInfo: SpanFitInformation,
        maxClippedLength: Int,
        maxClippedScore: Double
    ): Pair<Int, Int> {
        // Try to change left boundary
        val maxStart = start + maxClippedLength
        var clippedStart = start
        var step = CLIP_STEPS.size - 1
        while (step >= 0 && clippedStart <= maxStart) {
            val newStart = clippedStart + CLIP_STEPS[step]
            if (newStart > maxStart) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clipScore = fitInfo.score(ChromosomeRange(start, newStart, chromosome))
            if (clipScore < maxClippedScore) {
                clippedStart = newStart
                step = min(step + 1, CLIP_STEPS.size - 1)
            } else {
                step -= 1
            }
        }
        // Try to change right boundary
        val minEnd = end - maxClippedLength
        var clippedEnd = end
        step = CLIP_STEPS.size - 1
        while (step >= 0 && clippedEnd >= minEnd) {
            val newEnd = clippedEnd - CLIP_STEPS[step]
            if (newEnd < minEnd) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clipScore = fitInfo.score(ChromosomeRange(newEnd, end, chromosome))
            if (clipScore < maxClippedScore) {
                clippedEnd = newEnd
                step = min(step + 1, CLIP_STEPS.size - 1)
            } else {
                step -= 1
            }
        }
        return clippedStart to clippedEnd
    }


    /**
     * We want summary score to be robust wrt appending blocks of low significance,
     * so take into account top N% blocks p-values, otherwise we'll get fdr-blinking peaks, i.e.
     * peaks which are present for stronger fdr, but missing for more relaxed settings
     * Use length weighted mean to take into account difference in blocks lengths
     */
    private fun lengthWeightedScores(
        blocks: List<Range>,
        scores: List<Double>,
        fraction: Double = 0.5
    ): Double {
        require(blocks.size == scores.size) { "Different lengths of blocks and scores lists" }
        val sum = KahanSum()
        var l = 0
        blocks.zip(scores).sortedBy { it.second }
            .take(max(1, ceil(blocks.size * fraction).toInt()))
            .forEach { (b, p) ->
                sum += p * (b.startOffset - b.endOffset)
                l += b.startOffset - b.endOffset
            }
        return sum.result() / l
    }


    /**
     * Relaxed FDR should be useful for strict FDR only,
     * so limiting its application on relaxed FDR
     */
    fun relaxedLogFdr(logFdr: Double, relaxPower: Double = RELAX_POWER_DEFAULT) =
        max(logFdr, min(ln(0.1), logFdr * relaxPower))

    private const val RELAX_POWER_DEFAULT = 0.5

    private const val PSEUDO_COUNT: Int = 1

    private val CLIP_STEPS = intArrayOf(1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    private const val MAX_CLIPPED_LENGTH = 0.25
    private const val MAX_CLIPPED_DELTA = 0.5

    private val LOG_10 = ln(10.0)

    private val candidateBinsCounter = AtomicInteger()
    private val candidatePeaksCounter = AtomicInteger()
    private val resultPeaksCounter = AtomicInteger()

    private fun resetCounters() {
        candidateBinsCounter.set(0)
        candidatePeaksCounter.set(0)
        resultPeaksCounter.set(0)
    }

}