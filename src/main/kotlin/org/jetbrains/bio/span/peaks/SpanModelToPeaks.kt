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
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_MAX_LENGTH
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_CLIP_STEPS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_CLIP_MAX_SIGNAL
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FDR
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_HARD
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_LIGHT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_FRAGMENTATION_SPEED
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTATION_CHECKPOINT
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_FRAGMENTATION_MAX_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_MIN_SENSITIVITY
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SCORE_BLOCKS_GAP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SENSITIVITY_N
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SUMMITS_MIN_DISTANCE
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_SUMMITS_MIN_LENGTH
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


object SpanModelToPeaks {

    val LOG: Logger = LoggerFactory.getLogger(SpanModelToPeaks.javaClass)

    /**
     * Main method to compute peaks from the model.
     *
     * 1) Estimate posterior error probabilities PEP
     * 2) Find sensitivity setting as stable point wrt candidates number and average length.
     * 3) Compute required gap to compensate for extra high fragmentation, i.e. fraction of peaks when increasing gap
     * 4) Assign p-value to each peak based on combined p-values for cores (consequent foreground bins).
     *    In case when control track is present, we use Poisson CDF to estimate log P-value;
     *    otherwise, an average log PEP (posterior error probability) for bins in blocks is used.
     *    N% top significant blocks scores are aggregated using length-weighted average as P for peak.
     * 5) Compute qvalues by peaks p-values, filter by alpha.
     * 6) Fine-tuning boundaries of point-wise peaks according to the local signal.
     */
    fun getPeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        multipleTesting: MultipleTesting,
        sensitivityCmdArg: Double?,
        gapCmdArg: Int?,
        summits: Boolean,
        fragmentationLight: Double,
        fragmentationHard: Double,
        fragmentationSpeed: Double,
        clip: Double,
        blackListPath: Path? = null,
        parallel: Boolean = true,
        name: String? = null,
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
            sensitivityCmdArg to null
        } else {
            estimateSensitivity(
                genomeQuery, spanFitResults, logNullMembershipsMap, bitList2reuseMap, summits,
                parallel, name, cancellableState
            )
        }

        val blackList =
            if (blackListPath != null) {
                LOG.info("Loading blacklisted regions: $blackListPath")
                LocationsMergingList.load(genomeQuery, blackListPath)
            } else null

        LOG.info("Candidates selection with sensitivity: $sensitivity2use")

        val gap2use = when {
            gapCmdArg != null -> gapCmdArg
            summits -> 0
            else -> {
                LOG.info("${name ?: ""} Analysing fragmentation...")
                val candidateGapNs = IntArray(SPAN_FRAGMENTATION_MAX_GAP) {
                    estimateCandidatesNumberLens(
                        genomeQuery, fitInfo, logNullMembershipsMap, bitList2reuseMap,
                        sensitivity2use.first, it
                    ).n
                }
                estimateGap(candidateGapNs, name, fragmentationLight, fragmentationHard, fragmentationSpeed)
            }
        }
        LOG.info("${name ?: ""} Candidates selection with gap: $gap2use")

        val candidatesMap = genomeMap(genomeQuery, parallel = parallel) { chromosome ->
            cancellableState?.checkCanceled()
            if (!fitInfo.containsChromosomeInfo(chromosome)) {
                return@genomeMap emptyList<Range>()
            }
            val logNullMemberships = logNullMembershipsMap[chromosome]
            val bitList2reuse = bitList2reuseMap[chromosome]
            getChromosomeCandidates(
                chromosome, logNullMemberships, bitList2reuse,
                sensitivity2use.first, sensitivity2use.second, gap2use
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
        // Return empty list when nothing found
        if (candidatesTotal == 0) {
            return SpanPeaksResult(
                fdr, sensitivity2use.first, sensitivity2use.second, gap2use,
                genomeMap(genomeQuery) { emptyList() })
        }

        // Estimate signal and noise average signal by candidates
        val canEstimateSignalToNoise = clip > 0 && fitInfo is SpanAnalyzeFitInformation &&
                fitInfo.normalizedCoverageQueries?.all { it.areCachesPresent() } ?: false

        val (avgSignalDensity, avgNoiseDensity) = if (canEstimateSignalToNoise)
            estimateGenomeSignalNoiseAverage(genomeQuery, fitInfo, candidatesMap, parallel).apply {
                LOG.debug("${name ?: ""} Signal density $first, noise density $second")
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
        LOG.info("${name ?: ""} Adjusting pvalues ${multipleTesting.description}, N=${genomeLogPVals.length}")
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
                avgSignalDensity, avgNoiseDensity, clip,
                cancellableState = cancellableState
            )
        }
        return SpanPeaksResult(fdr, sensitivity2use.first, sensitivity2use.second, gap2use, peaks)
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
        candidatesNs: IntArray,
        name: String?,
        fragmentationLight: Double = SPAN_DEFAULT_FRAGMENTATION_LIGHT,
        fragmentationHard: Double = SPAN_DEFAULT_FRAGMENTATION_HARD,
        fragmentationSpeed: Double = SPAN_DEFAULT_FRAGMENTATION_SPEED,
    ): Int {
        val fragmentations = DoubleArray(candidatesNs.size) {
            candidatesNs[it].toDouble() / candidatesNs[0]
        }
        val tLightFragmentation = (1 until SPAN_FRAGMENTATION_CHECKPOINT).firstOrNull {
            fragmentations[it] <= fragmentationLight
        }
        LOG.debug(
            "${name ?: ""} Fragmentation @${SPAN_FRAGMENTATION_CHECKPOINT} = " +
                    "${fragmentations[SPAN_FRAGMENTATION_CHECKPOINT]}"
        )
        if (tLightFragmentation == null) {
            LOG.info("${name ?: ""} No fragmentation detected!")
            return 0
        } else {
            LOG.debug("${name ?: ""} Fragmentation light gap: $tLightFragmentation")
        }
        val speed = DoubleArray(fragmentations.size - 1) {
            fragmentations[it + 1] - fragmentations[it]
        }
        val tHardFragmentation = fragmentations.indices.firstOrNull {
            fragmentations[it] <= fragmentationHard
        }
        LOG.debug("${name ?: ""} Fragmentation hard gap: $tHardFragmentation")
        val tSpeedFragmentation = (tLightFragmentation / 2 until speed.size).firstOrNull {
            abs(speed[it]) < fragmentationSpeed
        }
        LOG.debug("${name ?: ""} Fragmentation speed gap: $tSpeedFragmentation")
        val finalFragmentationGap = when {
            tHardFragmentation != null && tSpeedFragmentation != null ->
                min(tHardFragmentation, tSpeedFragmentation)

            tHardFragmentation != null -> tHardFragmentation
            tSpeedFragmentation != null -> tSpeedFragmentation
            else -> tLightFragmentation
        }
        LOG.info("${name ?: ""} Fragmentation compensation gap: $finalFragmentationGap")
        return finalFragmentationGap
    }


    fun estimateSensitivity(
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>,
        summits: Boolean,
        parallel: Boolean,
        name: String?,
        cancellableState: CancellableState?
    ): Pair<Double, Double?> {
        cancellableState?.checkCanceled()

        LOG.info("${name ?: ""} Adjusting sensitivity...")
        val minLogNull = genomeQuery.get().minOf { logNullMembershipsMap[it].min() }
        // Limit value due to floating point errors
        val maxLogNull = min(
            SPAN_MIN_SENSITIVITY, genomeQuery.get().maxOf { logNullMembershipsMap[it].max() }
        )
        val sensitivities = linSpace(minLogNull, maxLogNull, SPAN_SENSITIVITY_N)
        // Compute candidates characteristics
        val (candidatesLogNs, candidatesLogALs) = getCandidatesCharacteristics(
            sensitivities,
            genomeQuery,
            spanFitResults,
            logNullMembershipsMap,
            bitList2reuseMap
        )
        val st = detectSensitivityTriangle(sensitivities, candidatesLogNs, candidatesLogALs)
        cancellableState?.checkCanceled()
        if (st != null) {
            LOG.info("${name ?: ""} Analysing candidates additive numbers...")
            val sensitivitiesLimited =
                sensitivities.slice(st.beforeMerge until st.stable).toDoubleArray()
            val (totals, news) = analyzeAdditiveCandidates(
                genomeQuery,
                spanFitResults.fitInfo,
                logNullMembershipsMap,
                bitList2reuseMap,
                sensitivitiesLimited,
                parallel
            )
            val newCandidatesList = sensitivitiesLimited.indices
                .map { news[it].toDouble() / totals[it].toDouble() }
            val minAdditionalIdx = newCandidatesList.indices.minByOrNull { newCandidatesList[it] }!!
            val minAdditionalSensitivity = sensitivitiesLimited[minAdditionalIdx]
            LOG.info("${name ?: ""} Minimal additional ${st.beforeMerge + minAdditionalIdx}: $minAdditionalSensitivity")
            return if (summits) {
                minAdditionalSensitivity to sensitivities[st.beforeMerge]
            } else {
                minAdditionalSensitivity to null
            }
        } else {
            LOG.error("${name ?: ""} Failed to estimate sensitivity, using defaults.")
            return ln(SPAN_DEFAULT_FDR) to null
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
        sensitivities: DoubleArray,
        candidatesLogNs: DoubleArray,
        candidatesLogALs: DoubleArray
    ): SensitivityInfo? {
        LOG.debug("Compute sensitivity triangle...")
        val n = sensitivities.size
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

    internal fun getCandidatesCharacteristics(
        sensitivities: DoubleArray,
        genomeQuery: GenomeQuery,
        spanFitResults: SpanFitResults,
        logNullMembershipsMap: GenomeMap<F64Array>,
        bitList2reuseMap: GenomeMap<BitList>
    ): Pair<DoubleArray, DoubleArray> {
        val n = sensitivities.size
        val candidatesLogNs = DoubleArray(n)
        val candidatesLogALs = DoubleArray(n)
        // println("Sensitivity\tGap\tCandidatesN\tCandidatesAL")
        for ((i, s) in sensitivities.withIndex()) {
            val ci = estimateCandidatesNumberLens(
                genomeQuery, spanFitResults.fitInfo, logNullMembershipsMap, bitList2reuseMap,
                s, 0
            )
            candidatesLogNs[i] = ln1p(ci.n.toDouble())
            candidatesLogALs[i] = ln1p(ci.averageLen)
            // println("$s\t0\t${ci.n}\t${ci.averageLen}")
        }
        return Pair(candidatesLogNs, candidatesLogALs)
    }

    internal fun linSpace(min: Double, max: Double, n: Int): DoubleArray {
        // If we use linear space here, we can't see the plots with merge
        // return DoubleArray(n) {
        //     min + (max - min) * it.toDouble() / (n - 1)
        // }
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
            getChromosomeCandidates(
                chromosome, logNullMemberships, bitList2reuse, sensitivities[0], null, 0,
            ).toRangeMergingList()
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
                getChromosomeCandidates(
                    chromosome, logNullMemberships, bitList2reuse, s, null, 0,
                ).toRangeMergingList()
            }
            val total = genomeQuery.get().sumOf { candidates[it].size }
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
                    chromosome, logNullMemberships, bitList2reuse, sensitivity, null, gap,
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
        sensitivitySummits: Double?,
        gap: Int
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
        val candidates = bins.aggregate(gap)
        if (sensitivitySummits == null) {
            return candidates
        }
        // Prepare to detect summits
        bins.set(0, bins.size(), false)
        for (i in 0 until logNullMemberships.length) {
            bins.set(i, logNullMemberships[i] <= sensitivitySummits)
        }
        // Keep summits at least 3 bin
        val summits = bins.aggregate().filter { it.length() >= SPAN_SUMMITS_MIN_LENGTH }
        bins.set(0, bins.size(), false)
        summits.forEach { (s, e) -> bins.set(s, e, true) }
        // merge summits by min relative distance
        var lastEnd = -1
        var lastD = Int.MAX_VALUE
        for ((start, end) in summits) {
            val d = (end - start) * SPAN_SUMMITS_MIN_DISTANCE
            if (lastEnd != -1) {
                if (start - lastEnd <= min(lastD, d)) {
                    bins.set(lastEnd, start)
                }
            }
            lastEnd = end
            lastD = d
        }
        // Add candidates not covered by more than one summit
        candidates.forEach { (s, e) ->
            if ((s until e).count { bins[it] && (it == s || !bins[it - 1]) } <= 1) {
                bins.set(s, e, true)
            }
        }
        // Final summits
        return bins.aggregate()
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
        clip: Double,
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
                fitInfo.normalizedCoverageQueries?.all { it.areCachesPresent() } ?: false
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
                    chromosome, start, end, fitInfo, avgSignalDensity, avgNoiseDensity, clipSignal = clip,
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
                    fitInfo.normalizedCoverageQueries?.all { it.areCachesPresent() } ?: false
        val readsControlAvailable =
            fitInfo is SpanAnalyzeFitInformation &&
                    fitInfo.isControlAvailable() &&
                    fitInfo.normalizedCoverageQueries?.all { it.areCachesPresent() } ?: false

        // Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
        // Replacing calls NormalizedCoverageQuery#score and NormalizedCoverageQuery#controlScore
        val treatmentCovs = when {
            readsTreatmentAvailable && fitInfo is SpanAnalyzeFitInformation ->
                fitInfo.normalizedCoverageQueries?.map { it.treatmentReads.get() } ?: emptyList()

            else -> emptyList()
        }

        val controlCovs: List<Coverage>
        val controlScales: List<Double>
        when {
            readsControlAvailable && fitInfo is SpanAnalyzeFitInformation -> {
                controlCovs = fitInfo.normalizedCoverageQueries!!.map { it.controlReads!!.get() }
                controlScales = fitInfo.normalizedCoverageQueries!!.map { it.coveragesNormalizedInfo.controlScale }
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
                fitInfo.normalizedCoverageQueries?.all { it.areCachesPresent() } ?: false
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
        clipSignal: Double = SPAN_DEFAULT_CLIP_MAX_SIGNAL,
        clipLength: Double = SPAN_CLIP_MAX_LENGTH,
        clipSteps: DoubleArray = SPAN_CLIP_STEPS,
    ): Pair<Int, Int> {
        if (avgSignalDensity <= avgNoiseDensity) {
            return start to end
        }
        // Additionally, clip peaks by local coverage signal
        val maxClippedDensity = avgNoiseDensity + clipSignal * (avgSignalDensity - avgNoiseDensity)
        val maxClippedSideLength = (end - start) * clipLength / 2
        val bin = fitInfo.binSize

        // Try to change the left boundary
        val maxStart = start + maxClippedSideLength
        var clippedStart = start
        var step = clipSteps.size - 1
        while (step >= 0 && clippedStart <= maxStart) {
            val newStart = clippedStart + (clipSteps[step] * bin).toInt()
            if (newStart > maxStart) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clippedDensity = fitInfo.score(ChromosomeRange(start, newStart, chromosome)) / (newStart - start)
            if (clippedDensity < maxClippedDensity) {
                clippedStart = newStart
                step = min(step + 1, clipSteps.size - 1)
            } else {
                step -= 1
            }
        }
        // Try to change the right boundary
        val minEnd = end - maxClippedSideLength
        var clippedEnd = end
        step = clipSteps.size - 1
        while (step >= 0 && clippedEnd >= minEnd) {
            val newEnd = clippedEnd - (clipSteps[step] * bin).toInt()
            if (newEnd < minEnd) {
                step -= 1
                continue
            }
            // Clip while clipped part score is less than average density
            val clippedDensity = fitInfo.score(ChromosomeRange(newEnd, end, chromosome)) / (end - newEnd)
            if (clippedDensity < maxClippedDensity) {
                clippedEnd = newEnd
                step = min(step + 1, clipSteps.size - 1)
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
