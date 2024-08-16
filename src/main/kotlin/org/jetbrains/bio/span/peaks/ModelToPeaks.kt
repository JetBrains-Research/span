package org.jetbrains.bio.span.peaks

import com.google.common.util.concurrent.AtomicDouble
import org.apache.commons.math3.stat.StatUtils
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.apache.commons.math3.util.Precision
import org.jetbrains.bio.dataframe.BitList
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.containers.toRangeMergingList
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.span.SpanResultsAnalysis
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_CHECKPOINT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_MAX_SHIFT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_NARROW_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_BROAD_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_STEPS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_LENGTH_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SIGNAL_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTATION_GAP_CHECKPOINT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_GAP_FRAGMENTATION_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SENSITIVITY_N
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.statistics.util.PoissonUtil
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.await
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import java.util.*
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong
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
        parallel: Boolean = true,
        cancellableState: CancellableState? = null,
    ): SpanPeaksResult {
        val fitInfo = spanFitResults.fitInfo
        // Prepare fit information for scores computations
        fitInfo.prepareData()

        // Collect candidates from model
        val logNullMembershipsMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                // HACK, it shouldn't be used
                // Cannot be null because of GenomeMap and empty F64Array is not supported
                return@genomeMap F64Array.full(1, 0.0)
            }
            getLogNulls(spanFitResults, chromosome)
        }
        val bitList2reuseMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            BitList(logNullMembershipsMap[chromosome].length)
        }

        val (sensitivity2use, gap2use) = estimateParams(
            genomeQuery, spanFitResults, logNullMembershipsMap, bitList2reuseMap,
            sensitivityCmdArg, gapCmdArg,
            blackListPath, parallel, cancellableState
        )

        val candidatesMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(
                chromosome, logNullMemberships, bitList2reuse, sensitivity2use, gap2use
            )
        }

        // Estimate signal and noise average signal by candidates
        val canEstimateSignalToNoise = clip && fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }

        val (avgSignalDensity, avgNoiseDensity) = if (canEstimateSignalToNoise)
            estimateGenomeSignalNoiseAverage(genomeQuery, fitInfo, candidatesMap, parallel).apply {
                LOG.debug("Signal density $first, noise density $second")
            }
        else
            null to null

        val readsControlAvailable = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.isControlAvailable() &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }

        // Collect peaks
        val peaks = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Peak>()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val candidates = candidatesMap[chromosome]
            val offsets = fitInfo.offsets(chromosome)
            getChromosomePeaksFromCandidates(
                chromosome, candidates, fitInfo, logNullMemberships, offsets,
                fdr, avgSignalDensity, avgNoiseDensity, readsControlAvailable,
                cancellableState = cancellableState
            )
        }
        return SpanPeaksResult(fdr, sensitivity2use, gap2use, peaks)
    }


    private fun estimateParams(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        blackListPath: Path?,
        parallel: Boolean,
        cancellableState: CancellableState?
    ): Pair<Double, Int> {
        if (sensitivityCmdArg != null && gapCmdArg != null) {
            return sensitivityCmdArg to gapCmdArg
        }
        cancellableState?.checkCanceled()

        LOG.info("Analysing sensitivity and gap...")
        val minLogNull = genomeQuery.get().minOf { logNullMembershipsMap[it].min() }
        // Limit value due to floating point errors
        val maxLogNull = min(-1e-10, genomeQuery.get().maxOf { logNullMembershipsMap[it].max() })
        val sensitivities = linSpace(minLogNull, maxLogNull)
        val si = detectSensitivityTriangle(
            genomeQuery, spanFitResults.fitInfo, logNullMembershipsMap, bitList2reuseMap, sensitivities
        )
        cancellableState?.checkCanceled()

        val fitInfo = spanFitResults.fitInfo

        LOG.info("Analysing autocorrelation...")
        val logNullPvals = getLogNullPvals(genomeQuery, spanFitResults, blackListPath)
        val autocorrelationScore = getAutocorrelation(logNullPvals, SPAN_AUTOCORRELATION_CHECKPOINT)
        LOG.info("Autocorrelation score: $autocorrelationScore")

        LOG.info("Analysing fragmentation...")
        val candidatesNGapI = estimateCandidatesNumberAvgLen(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
            sensitivities[si.stable], SPAN_FRAGMENTATION_GAP_CHECKPOINT
        ).first
        val candidatesNGap0 = estimateCandidatesNumberAvgLen(
            genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
            sensitivities[si.stable], 0
        ).first
        val fragmentationScore =
            if (candidatesNGap0 == 0) 1.0 else candidatesNGapI.toDouble() / candidatesNGap0
        LOG.info("Fragmentation score: $fragmentationScore")

        val (sensitivity2use, gap2use) = actualSensitivityGap(
            si, sensitivities, sensitivityCmdArg, gapCmdArg, autocorrelationScore, fragmentationScore,
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

    data class SensitivityInfo(
        val beforeMerge: Int, val stable: Int, val beforeNoise: Int,
    )

    fun detectSensitivityTriangle(
        genomeQuery: GenomeQuery,
        spanFitInformation: SpanFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivities: DoubleArray,
    ): SensitivityInfo {
        LOG.debug("Compute sensitivity triangle...")
        val n = sensitivities.size
        val candidatesNs = DoubleArray(n)
        val candidatesALs = DoubleArray(n)
        for ((i, s) in sensitivities.withIndex()) {
            val (cN, cAL) = estimateCandidatesNumberAvgLen(
                genomeQuery, spanFitInformation, logNullMembershipsMap, bitList2reuseMap, s, 0
            )
            candidatesNs[i] = ln1p(cN.toDouble())
            candidatesALs[i] = ln1p(cAL)
        }
        var maxArea = 0.0
        var i1 = -1
        var i2 = -1
        var i3 = -1
        for (i in 1 until n - 1) {
            val i1mab = findSensitivityTriangleMaxAreaBetween(
                candidatesNs, candidatesALs, 0, i, -1
            )
            val i3mab = findSensitivityTriangleMaxAreaBetween(
                candidatesNs, candidatesALs, i, n - 1, -1
            )
            if (i1mab.first == -1 || i3mab.first == -1) {
                continue
            }
            // We want both parts to be balanced so geometric mean optimization is better here
            val area = sqrt(i1mab.second * i3mab.second)
            if (area > maxArea) {
                maxArea = area
                i1 = i1mab.first
                i3 = i3mab.first
                i2 = i
            }
        }
        if (i1 == -1 || i2 == -1 || i3 == -1) {
            LOG.debug("Result beforeMerge: {}, stable: {}, beforeNoise: {}", i3, i2, i1)
            LOG.warn("Failed to estimate sensitivity triangle")
            return SensitivityInfo(0, (n - 1) / 2, n - 1)
        }
        // Update i3, i1 points to be closer to i2 for more accurate pivot estimations
        val i3mab = findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, i3, i2, -1)
        if (i3mab.first != -1) {
            i3 = i3mab.first
        }
        val i1mab = findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, i2, i1, -1)
        if (i1mab.first != -1) {
            i1 = i1mab.first
        }
        val result = SensitivityInfo(i1, i2, i3)
        LOG.debug(
            "Result beforeMerge: {}: {}, stable: {}: {}, beforeNoise: {}: {}",
            i1, sensitivities[i1], i2, sensitivities[i2], i3, sensitivities[i3],
        )
        return result
    }

    internal fun linSpace(min: Double, max: Double, n: Int = SPAN_SENSITIVITY_N): DoubleArray {
        // If we use linear space here, we can't see the plots with merge
//        return DoubleArray(n) {
//            min + (max - min) * it.toDouble() / (n - 1)
//        }
        require(min * max >= 0) { "Both min and max should have same sign, got $min, $max" }
        val sign = if (min + max < 0) -1 else 1
        val maxLog = log10(max * sign)
        val minLog = log10(min * sign)
        return DoubleArray(n) {
            sign * exp((minLog + (maxLog - minLog) * it.toDouble() / (n - 1)) * ln(10.0))
        }
    }

    private fun findSensitivityTriangleMaxAreaBetween(
        candidatesNs: DoubleArray,
        candidatesALs: DoubleArray,
        start: Int, end: Int, sign: Int = 1
    ): Pair<Int, Double> {
        if (start > end) {
            return findSensitivityTriangleMaxAreaBetween(candidatesNs, candidatesALs, end, start, sign)
        }
        var maxI = -1
        var maxArea = 0.0

        val startN = candidatesNs[start]
        val startAL = candidatesALs[start]
        val endN = candidatesNs[end]
        val endAL = candidatesALs[end]

        for (i in start + 1 until end) {
            val n = candidatesNs[i]
            val al = candidatesALs[i]
            var area = triangleSignedSquare(startN, startAL, n, al, endN, endAL)
            if (area * sign > 0)
                continue
            area = abs(area)
            if (area > maxArea) {
                maxI = i
                maxArea = area
            }
        }
        return maxI to maxArea
    }

    fun actualSensitivityGap(
        sensitivityInfos: SensitivityInfo,
        sensitivities: DoubleArray,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        autocorrelationScore: Double,
        fragmentationScore: Double,
    ): Pair<Double, Int> {
        val sensitivity2use = when {
            sensitivityCmdArg != null -> sensitivityCmdArg
            // Keep narrow signal narrow
            autocorrelationScore < SPAN_AUTOCORRELATION_NARROW_THRESHOLD ->
                sensitivities[sensitivityInfos.beforeMerge]
            // Stable point of HMM
            else -> sensitivities[sensitivityInfos.stable]
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
        spanFitInformation: SpanFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivities: DoubleArray,
        parallel: Boolean
    ): Pair<IntArray, IntArray> {
        // Collect candidates from model for most strict sensitivity
        // We assume that sensitivities are sorted ascending
        var candidatesPrev = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            if (!spanFitInformation.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>().toRangeMergingList()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(chromosome, logNullMemberships, bitList2reuse, sensitivities[0], 0)
                .toRangeMergingList()
        }

        val totals = IntArray(sensitivities.size)
        val news = IntArray(sensitivities.size)
        val candidatesPrevSize = genomeQuery.get().sumOf { candidatesPrev[it].size }
        totals[0] = candidatesPrevSize
        news[0] = candidatesPrevSize
        for (i in 1 until sensitivities.size) {
            val s = sensitivities[i]
            val candidates = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
                if (!spanFitInformation.containsChromosomeInfo(chromosome)) {
                    return@genomeMap emptyList<Range>().toRangeMergingList()
                }
                val logNullMemberships = logNullMembershipsMap[chromosome]
                val bitList2reuse = bitList2reuseMap[chromosome]
                getChromosomeCandidates(chromosome, logNullMemberships, bitList2reuse, s, 0).toRangeMergingList()
            }
            val candidatesSize = genomeQuery.get().sumOf { candidates[it].size }
            val total = candidatesSize
            val old = genomeQuery.get().sumOf { chromosome ->
                val chrCandidatesPrev = candidatesPrev[chromosome]
                candidates[chromosome].count { chrCandidatesPrev.intersectionLength(it) > 0 }
            }
            val new = total - old
            totals[i] = total
            news[i] = new
            candidatesPrev = candidates
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

    fun estimateCandidatesNumberAvgLen(
        genomeQuery: GenomeQuery,
        spanFitInformation: SpanFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        gap: Int
    ): Pair<Int, Double> {
        val candidatesN = AtomicInteger()
        val candidatesL = AtomicLong()
        genomeQuery.get().map { chromosome ->
            Callable {
                if (!spanFitInformation.containsChromosomeInfo(chromosome)) {
                    return@Callable
                }
                val logNullMemberships = logNullMembershipsMap[chromosome]
                val bitList2reuse = bitList2reuseMap[chromosome]
                val candidates = getChromosomeCandidates(
                    chromosome, logNullMemberships, bitList2reuse, sensitivity, gap
                )
                candidatesN.addAndGet(candidates.size)
                candidatesL.addAndGet(candidates.sumOf { it.length().toLong() })
            }
        }.await(true)
        return candidatesN.get() to
                (if (candidatesN.get() == 0) 0.0 else candidatesL.get().toDouble() / candidatesN.get())
    }


    fun getChromosomeCandidates(
        chromosome: Chromosome,
        logNullMemberships: F64Array,
        bitList2reuse: BitList?,
        sensitivity: Double,
        gap: Int,
    ): List<Range> {
        if ('_' in chromosome.name ||
            "random" in chromosome.name.lowercase() ||
            "un" in chromosome.name.lowercase()
        ) {
            LOG.trace("Ignore ${chromosome.name}: chromosome name looks like contig")
            return emptyList()
        }
        check(bitList2reuse == null || logNullMemberships.length == bitList2reuse.size())
        val bins = if (bitList2reuse != null)
            bitList2reuse
        else
            BitList(logNullMemberships.length)
        for (i in 0 until logNullMemberships.length) {
            bins.set(i, logNullMemberships[i] <= sensitivity)
        }
        return bins.aggregate(gap)
    }

    fun getLogNulls(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome
    ) = spanFitResults.logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)

    fun getChromosomePeaksFromCandidates(
        chromosome: Chromosome,
        candidates: List<Range>,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        fdr: Double,
        avgSignalDensity: Double?,
        avgNoiseDensity: Double?,
        readsControlAvailable: Boolean,
        cancellableState: CancellableState?
    ): List<Peak> {
        // Compute candidate bins and peaks with relaxed background settings
        // 1) Return broad peaks in case of broad modifications even for strict FDR settings
        // 2) Mitigate the problem when the number of peaks for strict FDR is much bigger than for relaxed FDR
        if (candidates.isEmpty()) {
            return emptyList()
        }

        // We want two invariants from peaks pvalues:
        // 1) The stricter FDR, the fewer peaks with smaller average length
        // 2) Peaks should not disappear when relaxing FDR
        // Peak score is computed as length-weighted average p-value in its consequent enriched bins.
        val peaksLogPvalues = estimateCandidatesLogPs(
            chromosome,
            candidates,
            fitInfo,
            logNullMemberships,
            offsets,
            readsControlAvailable,
            cancellableState
        )

        // Filter result peaks by Q values
        val peaksLogQValues = Fdr.qvalidate(peaksLogPvalues, logResults = true)
        val lnFdr = ln(fdr)
        val canEstimateScore = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        val resultPeaks = candidates.mapIndexedNotNull { i, (from, to) ->
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
        candidates: List<Range>,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        readsControlAvailable: Boolean,
        cancellableState: CancellableState?
    ): F64Array {
        // Compute significant blocks within candidate
        val candidatesBlocks = candidatesBlocks(candidates, logNullMemberships)
        // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
        // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
        val treatmentCovs = (fitInfo as SpanAnalyzeFitInformation)
            .normalizedCoverageQueries!!.map { it.treatmentReads.get() }
        val controlCovs: List<Coverage>
        val controlScales: List<Double>
        if (readsControlAvailable) {
            controlCovs = fitInfo.normalizedCoverageQueries!!.map { it.controlReads!!.get() }
            controlScales = fitInfo.normalizedCoverageQueries!!
                .map { it.coveragesNormalizedInfo.controlScale }
        } else {
            controlCovs = emptyList()
            controlScales = emptyList()
        }

        val peaksLogPvalues = F64Array(candidates.size) { idx ->
            cancellableState?.checkCanceled()
            var blocks = candidatesBlocks[idx]
            // Use the whole peaks for estimation may produce more significant results
            // beneficial for short peaks, TFs, ATAC-seq etc.
            if (blocks.size <= 1) {
                blocks = listOf(candidates[idx])
            }
            val blocksLogPs = blocks.map { (from, to) ->
                // Model posterior log error probability for block
                val modelLogPs = (from until to).sumOf { logNullMemberships[it] }
                if (!readsControlAvailable) {
                    return@map modelLogPs
                }
                val blockStart = offsets[from]
                val blockEnd = if (to < offsets.size) offsets[to] else chromosome.length
                val chromosomeRange = ChromosomeRange(blockStart, blockEnd, chromosome)
                // val score = fitInfo.score(chromosomeRange)
                val score = SpanAnalyzeFitInformation.fastScore(treatmentCovs, chromosomeRange)
                // val controlScore = fitInfo.controlScore(chromosomeRange)
                val controlScore =
                    SpanAnalyzeFitInformation.fastControlScore(controlCovs, controlScales, chromosomeRange)
                // Use +1 as a pseudo count to compute Poisson CDF
                val signalLogPs = PoissonUtil.logPoissonCdf(
                    ceil(score).toInt() + 1,
                    controlScore + 1
                )
                // Combine both model and signal estimations
                return@map (modelLogPs + signalLogPs) / 2
            }
            return@F64Array lengthWeightedScores(blocks, blocksLogPs)
        }
        return peaksLogPvalues
    }

    fun candidatesBlocks(
        candidates: List<Range>,
        logNullMemberships: F64Array
    ) = candidates.map { (start, end) ->
        val p = StatUtils.percentile(
            logNullMemberships.slice(start, end).toDoubleArray(), SPAN_SCORE_BLOCKS * 100
        )
        return@map BitList(end - start) {
            logNullMemberships[it + start] <= p
        }.aggregate(SPAN_DEFAULT_GAP).map { Range(start + it.startOffset, start + it.endOffset) }
    }

    fun estimateGenomeSignalNoiseAverage(
        genomeQuery: GenomeQuery,
        fitInfo: SpanFitInformation,
        candidates: GenomeMap<List<Range>>,
        parallel: Boolean
    ): Pair<Double, Double> {
        val canEstimateSignalToNoise = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        if (!canEstimateSignalToNoise) {
            return 0.0 to 0.0
        }
        val sumSignalScoreA = AtomicDouble()
        val sumSignalLengthA = AtomicLong()
        val sumNoiseScoreA = AtomicDouble()
        val sumNoiseLengthA = AtomicLong()
        // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
        // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
        val treatmentCovs = (fitInfo as SpanAnalyzeFitInformation)
            .normalizedCoverageQueries!!.map { it.treatmentReads.get() }
        genomeQuery.get().map {  chromosome ->
            Callable {
                if (!fitInfo.containsChromosomeInfo(chromosome) || chromosome !in candidates) {
                    return@Callable
                }
                val chrCandidates = candidates[chromosome]
                val offsets = fitInfo.offsets(chromosome)
                var prevNoiseStart = 0
                chrCandidates.forEach { (from, to) ->
                    val start = offsets[from]
                    val end = if (to < offsets.size) offsets[to] else chromosome.length
                    val range = ChromosomeRange(start, end, chromosome)
                    // val score = fitInfo.score(range)
                    val score = SpanAnalyzeFitInformation.fastScore(treatmentCovs, range)
                    sumSignalScoreA.addAndGet(score)
                    sumSignalLengthA.addAndGet(end.toLong() - start)
                    val rangeNoise = ChromosomeRange(prevNoiseStart, start, chromosome)
                    // val scoreNoise = fitInfo.score(rangeNoise)
                    val scoreNoise = SpanAnalyzeFitInformation.fastScore(treatmentCovs, rangeNoise)
                    sumNoiseScoreA.addAndGet(scoreNoise)
                    sumNoiseLengthA.addAndGet(start.toLong() - prevNoiseStart)
                    prevNoiseStart = end
                }
                if (prevNoiseStart < chromosome.length) {
                    val rangeNoise = ChromosomeRange(prevNoiseStart, chromosome.length, chromosome)
                    // val scoreNoise = fitInfo.score(rangeNoise)
                    val scoreNoise = SpanAnalyzeFitInformation.fastScore(treatmentCovs, rangeNoise)
                    sumNoiseScoreA.addAndGet(scoreNoise)
                    sumNoiseLengthA.addAndGet(chromosome.length.toLong() - prevNoiseStart)
                }
            }
        }.await(parallel)
        val sumSignalScore = sumSignalScoreA.get()
        val sumSignalLength = sumSignalLengthA.get()
        val sumNoiseScore = sumNoiseScoreA.get()
        val sumNoiseLength = sumNoiseLengthA.get()
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