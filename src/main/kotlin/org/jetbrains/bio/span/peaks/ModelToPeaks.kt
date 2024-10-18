package org.jetbrains.bio.span.peaks

import com.google.common.util.concurrent.AtomicDouble
import gnu.trove.list.array.TDoubleArrayList
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
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_AUTOCORRELATION_MAX_SHIFT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_BROAD_AC_MIN_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_BROAD_EXTRA_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_STEPS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTATION_MAX_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTED_EXTRA_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTED_MAX_THRESHOLD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_LENGTH_CLIP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MIN_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SENSITIVITY_N
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SIGNAL_CLIP
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
import java.util.concurrent.atomic.AtomicLong
import kotlin.math.*


object ModelToPeaks {

    val LOG: Logger = LoggerFactory.getLogger(ModelToPeaks.javaClass)

    /**
     * Main method to compute peaks from the model.
     *
     * 1) Estimate posterior probabilities.
     * 2) Find sensitivity setting as stable point wrt candidates number and average length.
     * 3) Assign p-value to each peak based on combined p-values for cores (consequent foreground bins).
     *    In case when control track is present, we use Poisson CDF to estimate log P-value;
     *    otherwise, an average log PEP (posterior error probability) for bins in blocks is used.
     *    N% top significant blocks scores are aggregated using length-weighted average as P for peak.
     * 4) Compute qvalues by peaks p-values, filter by alpha.
     * 5) Optional clipping to fine-tune boundaries of point-wise peaks according to the local signal.
     */
    fun getPeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        multipleTesting: MultipleTesting,
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
                return@genomeMap FAKE_ARRAY
            }
            getLogNulls(spanFitResults, chromosome)
        }
        val bitList2reuseMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            BitList(logNullMembershipsMap[chromosome].length)
        }
        val sensitivity2use = if (sensitivityCmdArg != null) {
            sensitivityCmdArg
        } else {
            // We want to be able to get shorter peaks with stringent fdr values
            min(ln(fdr), estimateSensitivity(
                genomeQuery, spanFitResults, logNullMembershipsMap, bitList2reuseMap,
                parallel, cancellableState
            ))
        }

        val blackList =
            if (blackListPath != null) {
                LOG.info("Loading blacklisted regions: $blackListPath")
                LocationsMergingList.load(genomeQuery, blackListPath)
            } else null

        LOG.info("Candidates selection with sensitivity: $sensitivity2use")

        val gap2use = if (gapCmdArg != null) {
            gapCmdArg
        } else {
            LOG.info("Analysing pvalues autocorrelation...")
            val logNullPvals = getLogNullPvals(genomeQuery, spanFitResults, blackList)
            val logNullPValsCorrelations = computeCorrelations(logNullPvals)
            val avgAutoCorrelation = logNullPValsCorrelations.average()
            LOG.info("Average autocorrelation score: $avgAutoCorrelation")

            LOG.info("Analysing fragmentation...")
            val candidateGapNs = IntArray(SPAN_FRAGMENTATION_MAX_GAP) {
                estimateCandidatesNumberLens(
                    genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                    sensitivity2use, it
                ).n
            }
            val avgFragmentationScore = candidateGapNs.sumOf { it.toDouble() / candidateGapNs[0] }
            LOG.info("Average fragmentation score: $avgFragmentationScore")

            estimateGap(avgAutoCorrelation, avgFragmentationScore)
        }
        LOG.info("Candidates selection with gap: $gap2use")

        val candidatesMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(
                chromosome, logNullMemberships, bitList2reuse,
                sensitivity2use, gap2use
            )
        }

        // Estimate total number of candidates and compute indices by chromosome
        var candidatesTotal = 0
        val chromosomeCandidatesOffsets = hashMapOf<Chromosome, Pair<Int, Int>>()
        genomeQuery.get().forEach { chromosome ->
            val chrCandidatesN = candidatesMap[chromosome].size
            chromosomeCandidatesOffsets[chromosome] = candidatesTotal to candidatesTotal + chrCandidatesN
            candidatesTotal += chrCandidatesN
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

        val logPValsMap = collectPVals(
            genomeQuery, fitInfo, candidatesMap, logNullMembershipsMap,
            parallel, cancellableState
        )

        // Adjust pvalues globally
        val genomeLogPVals = F64Array(candidatesTotal)
        genomeQuery.get().forEach { chromosome ->
            val (start, end) = chromosomeCandidatesOffsets[chromosome]!!
            if (start == end) {
                return@forEach
            }
            val chrLogPVals = logPValsMap[chromosome]
            check(chrLogPVals.length == end - start)
            for (i in start until end) {
                genomeLogPVals[i] = chrLogPVals[i - start]
            }
        }

        // Adjust globally log pValues -> log qValues
        LOG.info("Adjusting pvalues ${multipleTesting.description}, N=${genomeLogPVals.length}")
        val genomeLogQVals = if (multipleTesting == MultipleTesting.BH)
            Fdr.qvalidate(genomeLogPVals, logResults = true)
        else
            F64Array(genomeLogPVals.length) { genomeLogPVals[it] + ln(genomeLogPVals.length.toDouble()) }

        // Collect peaks together from all chromosomes
        val peaks = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Peak>()
            }
            val candidates = candidatesMap[chromosome]
            if (candidates.isEmpty()) {
                return@genomeMap emptyList<Peak>()
            }
            val offsets = fitInfo.offsets(chromosome)
            val (start, end) = chromosomeCandidatesOffsets[chromosome]!!
            val logPVals = genomeLogPVals.slice(start, end)
            val logQVals = genomeLogQVals.slice(start, end)
            getChromosomePeaksFromCandidates(
                chromosome, candidates, fitInfo, offsets,
                logPVals, logQVals,
                fdr, blackList,
                avgSignalDensity, avgNoiseDensity,
                cancellableState = cancellableState
            )
        }
        return SpanPeaksResult(fdr, sensitivity2use, gap2use, peaks)
    }

    private fun collectPVals(
        genomeQuery: GenomeQuery,
        fitInfo: SpanFitInformation,
        candidatesMap: GenomeMap<List<Range>>,
        logNullMembershipsMap: GenomeMap<F64Array>,
        parallel: Boolean,
        cancellableState: CancellableState?
    ): GenomeMap<F64Array> {
        val logPValsMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                // Return fake array to keep typing
                return@genomeMap FAKE_ARRAY
            }
            val candidates = candidatesMap[chromosome]
            if (candidates.isEmpty()) {
                return@genomeMap FAKE_ARRAY
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val offsets = fitInfo.offsets(chromosome)
            return@genomeMap estimateCandidatesLogPs(
                chromosome,
                candidates,
                fitInfo,
                logNullMemberships,
                offsets,
                cancellableState
            )
        }
        return logPValsMap
    }

    fun estimateGap(
        avgAutocorrelation: Double,
        avgFragmentation: Double
    ): Int {
        LOG.debug("Estimating gap...")
        var extraGap = 0
        if (avgAutocorrelation > SPAN_BROAD_AC_MIN_THRESHOLD) {
            LOG.info(
                "Autocorrelation $avgAutocorrelation > $SPAN_BROAD_AC_MIN_THRESHOLD, " +
                        "adding +$SPAN_BROAD_EXTRA_GAP"
            )
            extraGap += SPAN_BROAD_EXTRA_GAP
        }
        if (avgFragmentation < SPAN_FRAGMENTED_MAX_THRESHOLD) {
            LOG.info(
                "Fragmentation $avgFragmentation < $SPAN_FRAGMENTED_MAX_THRESHOLD, " +
                        "adding +$SPAN_FRAGMENTED_EXTRA_GAP"
            )
            extraGap += SPAN_FRAGMENTED_EXTRA_GAP
        }
        return SPAN_DEFAULT_GAP + extraGap
    }


    fun estimateSensitivity(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        parallel: Boolean,
        cancellableState: CancellableState?
    ): Double {
        cancellableState?.checkCanceled()

        LOG.info("Adjusting sensitivity...")
        val minLogNull = genomeQuery.get().minOf { logNullMembershipsMap[it].min() }
        // Limit value due to floating point errors
        val maxLogNull = min(SPAN_MIN_SENSITIVITY, genomeQuery.get().maxOf { logNullMembershipsMap[it].max() })
        val sensitivities = linSpace(minLogNull, maxLogNull, SPAN_SENSITIVITY_N)
        val si = detectSensitivityTriangle(
            genomeQuery, spanFitResults, logNullMembershipsMap, bitList2reuseMap, sensitivities
        )
        cancellableState?.checkCanceled()
        if (si != null) {
            LOG.info("Analysing candidates additive numbers...")
            val sensitivitiesLimited =
                sensitivities.slice(si.beforeMerge until si.stable).toDoubleArray()
            val (totals, news) = analyzeAdditiveCandidates(
                genomeQuery,
                spanFitResults.fitInfo,
                logNullMembershipsMap,
                bitList2reuseMap,
                sensitivitiesLimited,
                parallel
            )
            val minAdditionalIdx = sensitivitiesLimited.indices
                .minByOrNull { news[it].toDouble() / totals[it].toDouble() }!!
            val minAdditionalSensitivity = sensitivitiesLimited[minAdditionalIdx]
            LOG.info("Minimal additional ${si.beforeMerge + minAdditionalIdx}: $minAdditionalSensitivity")
            return minAdditionalSensitivity
        } else {
            LOG.error("Failed to estimate sensitivity, using defaults.")
            return SPAN_DEFAULT_SENSITIVITY
        }
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

    /**
     * Detects major semantic changes along sensitivity values.
     * 1) Before merge adjacent candidates
     * 2) Merge - stable
     * 3) Before noise calling
     * See [SensitivityInfo] for details
     */
    fun detectSensitivityTriangle(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivities: DoubleArray,
    ): SensitivityInfo? {
        LOG.debug("Compute sensitivity triangle...")
        val n = sensitivities.size
        val candidatesLogNs = DoubleArray(n)
        val candidatesLogALs = DoubleArray(n)
//        println("Sensitivity\tGap\tCandidatesN\tCandidatesAL")
        for ((i, s) in sensitivities.withIndex()) {
            val ci = estimateCandidatesNumberLens(
                genomeQuery, spanFitResults.fitInfo, logNullMembershipsMap, bitList2reuseMap,
                s, 0
            )
            candidatesLogNs[i] = ln1p(ci.n.toDouble())
            candidatesLogALs[i] = ln1p(ci.averageLen)
//            println("$s\t0\t${ci.n}\t${ci.averageLen}")
        }
        var maxArea = 0.0
        var i1 = -1
        var i2 = -1
        var i3 = -1
        for (i in (n * 0.2).toInt()..(n * 0.8).toInt()) {
            val i1mab = findSensitivityTriangleMaxAreaBetween(
                candidatesLogNs, candidatesLogALs, 0, i, -1
            )
            val i3mab = findSensitivityTriangleMaxAreaBetween(
                candidatesLogNs, candidatesLogALs, i, n - 1, -1
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
            return null
        }
        // Update i3, i1 points to be closer to i2 for more accurate pivot estimations
        val i3mab = findSensitivityTriangleMaxAreaBetween(candidatesLogNs, candidatesLogALs, i3, i2, -1)
        if (i3mab.first != -1) {
            i3 = i3mab.first
        }
        val i1mab = findSensitivityTriangleMaxAreaBetween(candidatesLogNs, candidatesLogALs, i2, i1, -1)
        if (i1mab.first != -1) {
            i1 = i1mab.first
        }
        val result = SensitivityInfo(i1, i2, i3)
        LOG.info(
            "Result beforeMerge: {}: {}, stable: {}: {}, beforeNoise: {}: {}",
            i1, sensitivities[i1], i2, sensitivities[i2], i3, sensitivities[i3],
        )
        return result
    }

    internal fun linSpace(min: Double, max: Double, n: Int): DoubleArray {
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
        blackList: LocationsMergingList?
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
            .sumOf { spanFitResults.logNullMemberships[it.name]!!.f64Array(SpanModelFitExperiment.NULL).length }
        val result = DoubleArray(totalBins) { 0.0 }
        var i = 0
        var blackListIgnored = 0
        for (chr in limitedQuery.get()) {
            val logNullPeps = spanFitResults.logNullMemberships[chr.name]!!.f64Array(SpanModelFitExperiment.NULL)
            for (j in 0 until logNullPeps.length) {
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
            LOG.debug("Marked {} / {} blacklisted regions", blackListIgnored, totalBins)
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

    private fun getAutocorrelation(values: DoubleArray, d: Int): Double {
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

    data class CandidatesInfo(
        val n: Int, val averageLen: Double, val medianLen: Double, val maxLen: Int
    )

    fun estimateCandidatesNumberLens(
        genomeQuery: GenomeQuery,
        spanFitInformation: SpanFitInformation,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        sensitivity: Double,
        gap: Int
    ): CandidatesInfo {
        val lens = TDoubleArrayList()
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
                if (candidates.isNotEmpty()) {
                    synchronized(lens) {
                        lens.add(DoubleArray(candidates.size) { candidates[it].length().toDouble() })
                    }
                }
            }
        }.await(true)
        return CandidatesInfo(
            lens.size(),
            if (lens.size() == 0) 0.0 else lens.sum() / lens.size(),
            if (lens.size() == 0) 0.0 else StatUtils.percentile(lens.toArray(), 50.0),
            if (lens.size() == 0) 0 else lens.max().toInt(),
        )
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
        offsets: IntArray,
        logPVals: F64Array,
        logQVals: F64Array,
        fdr: Double,
        blackList: LocationsMergingList?,
        avgSignalDensity: Double?,
        avgNoiseDensity: Double?,
        cancellableState: CancellableState?
    ): List<Peak> {
        // Compute candidate bins and peaks with relaxed background settings
        // 1) Return broad peaks in case of broad modifications even for strict FDR settings
        // 2) Mitigate the problem when the number of peaks for strict FDR is much bigger than for relaxed FDR
        if (candidates.isEmpty()) {
            return emptyList()
        }

        val lnFdr = ln(fdr)
        val canEstimateScore = fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        val resultPeaks = candidates.mapIndexedNotNull { i, (from, to) ->
            cancellableState?.checkCanceled()
            val logPValue = logPVals[i]
            val logQValue = logQVals[i]
            if (logPValue > lnFdr || logQValue > lnFdr) {
                return@mapIndexedNotNull null
            }
            var start = offsets[from]
            var end = if (to < offsets.size) offsets[to] else chromosome.length
            if (blackList != null && blackList.intersects(Location(start, end, chromosome))) {
                return@mapIndexedNotNull null
            }
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

    /**
     * We want two invariants from peaks pvalues:
     * 1) The stricter FDR, the fewer peaks with smaller average length
     * 2) Peaks should not disappear when relaxing FDR
     * Peak score is computed as length-weighted average p-value in its consequent enriched bins.
     */
    fun estimateCandidatesLogPs(
        chromosome: Chromosome,
        candidates: List<Range>,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
        cancellableState: CancellableState?
    ): F64Array {
        // TODO[oleg] support SpanCompareFitInformation
        val readsTreatmentAvailable =
            fitInfo is SpanAnalyzeFitInformation &&
                    fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }
        val readsControlAvailable =
            fitInfo is SpanAnalyzeFitInformation &&
                    fitInfo.isControlAvailable() &&
                    fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }

        // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
        // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
        val treatmentCovs = when {
            readsTreatmentAvailable && fitInfo is SpanAnalyzeFitInformation ->
                fitInfo.normalizedCoverageQueries!!.map { it.treatmentReads.get() }

            else -> emptyList()
        }

        val controlCovs: List<Coverage>
        val controlScales: List<Double>
        when {
            readsControlAvailable && fitInfo is SpanAnalyzeFitInformation -> {
                controlCovs = fitInfo.normalizedCoverageQueries!!.map { it.controlReads!!.get() }
                controlScales = fitInfo.normalizedCoverageQueries!!
                    .map { it.coveragesNormalizedInfo.controlScale }
            }

            else -> {
                controlCovs = emptyList()
                controlScales = emptyList()
            }
        }
        val peaksLogPvalues = F64Array(candidates.size) { idx ->
            cancellableState?.checkCanceled()
            val candidate = candidates[idx]
            // Compute significant blocks within candidate
            var blocks = candidateBlocks(logNullMemberships, candidate.startOffset, candidate.endOffset)
            // Expand single block to the whole candidate, useful for short peaks
            if (blocks.size == 1) {
                blocks = listOf(Range(candidate.startOffset, candidate.endOffset))
            }
            val blocksLogPs = blocks.map { (from, to) ->
                // Model posterior log error probability for block
                val modelLogPs = (from until to).sumOf { logNullMemberships[it] }
                if (!readsTreatmentAvailable || !readsControlAvailable) {
                    return@map modelLogPs
                }
                val blockStart = offsets[from]
                val blockEnd = if (to < offsets.size) offsets[to] else chromosome.length
                val chromosomeRange = ChromosomeRange(blockStart, blockEnd, chromosome)
                // val score = fitInfo.score(chromosomeRange)
                val score = SpanAnalyzeFitInformation.fastScore(treatmentCovs, chromosomeRange)
                // val controlScore = fitInfo.controlScore(chromosomeRange)
                val controlScore = SpanAnalyzeFitInformation.fastControlScore(
                    controlCovs, controlScales, chromosomeRange
                )
                // Use +1 as a pseudo count to compute Poisson CDF
                val signalLogPs = PoissonUtil.logPoissonCdf(
                    ceil(score).toInt() + 1,
                    controlScore + 1
                )
                // Combine both model and signal estimations
                return@map -sqrt(modelLogPs * signalLogPs)
            }
            return@F64Array lengthWeightedScores(blocks, blocksLogPs)
        }
        return peaksLogPvalues
    }

    fun candidateBlocks(
        logNullMemberships: F64Array,
        start: Int,
        end: Int
    ): List<Range> {
        val p = StatUtils.percentile(
            logNullMemberships.slice(start, end).toDoubleArray(), SPAN_SCORE_BLOCKS * 100
        )
        return BitList(end - start) {
            logNullMemberships[it + start] <= p
        }.aggregate(SPAN_SCORE_BLOCKS_GAP).map { Range(start + it.startOffset, start + it.endOffset) }
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
        genomeQuery.get().map { chromosome ->
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
        clipSignal: Double = SPAN_SIGNAL_CLIP,
        clipLength: Double = SPAN_LENGTH_CLIP,
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

    private val FAKE_ARRAY = F64Array.full(1, 0.0)

}