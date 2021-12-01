package org.jetbrains.bio.span

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.experiments.fit.SpanFitResults
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment
import org.jetbrains.bio.experiments.fit.f64Array
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.BenjaminiHochberg
import org.jetbrains.bio.statistics.hypothesis.StofferLiptakTest
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.log10
import kotlin.math.max
import kotlin.math.min

/**
 * The islands are called in five steps.
 *
 * 1) Estimate posterior probabilities
 * 2) Pick candidate bins with probability of H_0 <= alpha * [multiplier]
 * 3) Using gap merge bins into candidate islands
 * 4) Assign p-value to each island using Stoffer-Liptak test
 * 5) Compute qvalues on islands p-values, filter by alpha
 *
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 * @param multiplier is used to
 * 1) Return broad peaks in case of broad modifications even for strict FDR settings
 * 2) Mitigate the problem when number of peaks for strict FDR is much bigger than for relaxed FDR
 * @param noclip Do not clip islands to increase density when true
 */
fun SpanFitResults.getIslands(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    multiplier: Double = 1e2,
    noclip: Boolean = false,
    cancellableState: CancellableState? = null
): List<Peak> {
    val progress = Progress { title = "Computing islands" }.bounded(genomeQuery.get().size.toLong())
    fitInfo.prepareScores()
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        val chromosomePeaks = getChromosomePeaks(chromosome, fdr, gap, multiplier, noclip)
        progress.report(1)
        chromosomePeaks
    }
    progress.done()
    return genomeQuery.get().flatMap { map[it] }
}

private fun SpanFitResults.getChromosomePeaks(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    multiplier: Double,
    noclip: Boolean
): List<Peak> {
    // Check that we have information for requested chromosome
    val chromosomePeaks = if (chromosome.name in fitInfo.chromosomesSizes) {
        getChromosomeIslands(
            logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
            fitInfo.offsets(chromosome),
            chromosome,
            fdr,
            gap,
            multiplier,
            noclip
        )
    } else {
        SpanFitResults.LOG.debug(
            "NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}"
        )
        emptyList()
    }
    return chromosomePeaks
}

internal fun SpanFitResults.getChromosomeIslands(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    multiplier: Double,
    noclip: Boolean
): List<Peak> {
    // Compute candidate bins and islands
    val nullMemberships = logNullMemberships.exp()
    val candidateBins = BitterSet(logNullMemberships.size).apply {
        0.until(size()).filter { nullMemberships[it] <= min(0.1, fdr * multiplier) }.forEach(::set)
    }
    Peak.LOG.debug(
        "$chromosome: candidate bins ${candidateBins.cardinality()}/${logNullMemberships.size}"
    )
    val candidateIslands = candidateBins.aggregate(gap)
    val filteredIslands = candidateIslands.filter { (i, j) ->
        (i until j).any { nullMemberships[it] <= fdr }
    }
    if (filteredIslands.isEmpty()) {
        return emptyList()
    }
    val stofferLiptakTest = ISLANDS_STOFFER_LIPTAK_CACHE.get(this to chromosome) {
        StofferLiptakTest(nullMemberships.data)
    }
    // Apply Stoffer-Liptak test to correct dependence between consequent p-values
    val islandsPValues = F64Array(filteredIslands.size) { islandIndex ->
        val (i, j) = filteredIslands[islandIndex]
        val pValues = (i until j).map { nullMemberships[it] }.filter { it <= fdr }.toDoubleArray()
        stofferLiptakTest.combine(pValues)
    }

    val islandQValues = BenjaminiHochberg.adjust(islandsPValues)
    val resultIslands = filteredIslands.indices.filter { islandQValues[it] < fdr }
    var clipStart = 0L
    var clipEnd = 0L
    val clippedIslands = resultIslands.indices
        .map { islandIndex ->
            val (i, j) = filteredIslands[islandIndex]
            val start = offsets[i]
            val end = if (j < offsets.size) offsets[j] else chromosome.length
            // Optimize length
            val (clippedStart, clippedEnd) = if (noclip)
                start to end
            else {
                // Threshold limiting max clipped out density, more strict fdr should allow more clipping
                // Fdr      DensityMultiplier
                // 1e-1     0.2
                // 1e-2     0.6
                // 1e-3     0.73
                // 1e-4     0.8
                // 1e-6     0.86
                // 1e-10    0.92
                val densityMultiplier = 0.2 + 0.8 * max(0.0, 1 + 1 / log10(fdr))
                val maxClip = max(0, ((end - start) - (end - start) / (j - i) * gap / 2) / 2)
                clipIsland(start, end, maxClip, densityMultiplier) { s: Int, e: Int ->
                    fitInfo.score(ChromosomeRange(s, e, chromosome)) / (e - s)
                }
            }
            clipStart += (clippedStart - start)
            clipEnd += (end - clippedEnd)
            Peak(
                chromosome = chromosome,
                startOffset = clippedStart,
                endOffset = clippedEnd,
                mlogpvalue = -log10(islandsPValues[islandIndex]),
                mlogqvalue = -log10(islandQValues[islandIndex]),
                // Value is either coverage of fold change
                value = fitInfo.score(ChromosomeRange(clippedStart, clippedEnd, chromosome)),
                // Score should be proportional original q-value
                score = min(1000.0, -10 * log10(islandQValues[islandIndex])).toInt()
            )
        }
    Peak.LOG.debug(
        "$chromosome: islands result/filtered/candidate " +
                "${clippedIslands.size}/${filteredIslands.size}/${candidateIslands.size} " +
                "average clip start/end " +
                "${clipStart.toDouble() / max(1, clippedIslands.size)}/${
                    clipEnd.toDouble() / max(1, clippedIslands.size)
                }"
    )
    return clippedIslands
}


private val CLIP_STEPS = intArrayOf(1, 2, 3, 5, 10, 20, 50)

/**
 * Tries to reduce range by [CLIP_STEPS] from both sides while increasing [density].
 *
 * @param start Initial start of island
 * @param end Initial end of island
 * @param maxClip Limits maximum length to be clipped
 * @param multiplier Limits maximum density of clipped out fragment as initial island [density] * [multiplier]
 * @param density Density function, i.e. coverage / length
 */
internal fun clipIsland(
    start: Int,
    end: Int,
    maxClip: Int,
    multiplier: Double,
    density: (Int, Int) -> (Double)
): Pair<Int, Int> {
    val maxClippedDensity = density(start, end) * multiplier
    // Try to change left boundary
    val maxStart = start + maxClip
    var currentStart = start
    var step = CLIP_STEPS.size - 1
    while (step >= 0 && currentStart <= maxStart) {
        val newStart = currentStart + CLIP_STEPS[step]
        if (newStart > maxStart) {
            step -= 1
            continue
        }
        // Clip while clipped part density is less than average density
        val newDensity = density(start, newStart)
        if (newDensity < maxClippedDensity) {
            currentStart = newStart
            step = min(step + 1, CLIP_STEPS.size - 1)
        } else {
            step -= 1
        }
    }
    // Try to change right boundary
    val minEnd = end - maxClip
    var currentEnd = end
    step = CLIP_STEPS.size - 1
    while (step >= 0 && currentEnd >= minEnd) {
        val newEnd = currentEnd - CLIP_STEPS[step]
        if (newEnd < minEnd) {
            step -= 1
            continue
        }
        // Clip while clipped part density is less than average density
        val newDensity = density(newEnd, end)
        if (newDensity < maxClippedDensity) {
            currentEnd = newEnd
            step = min(step + 1, CLIP_STEPS.size - 1)
        } else {
            step -= 1
        }
    }
    return currentStart to currentEnd
}

/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * Using caches with weak values to avoid memory overflow.
 */
private val ISLANDS_STOFFER_LIPTAK_CACHE: Cache<Pair<SpanFitResults, Chromosome>, StofferLiptakTest> =
    CacheBuilder.newBuilder().weakValues().build()