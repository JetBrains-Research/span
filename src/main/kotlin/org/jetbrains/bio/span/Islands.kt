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
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.log10
import kotlin.math.min

/**
 * The islands are called in five steps.
 *
 * 1) Estimate posterior probabilities
 * 2) Pick candidate bins with probability of H_0 <= alpha * multiplier
 * 3) Using gap merge bins into candidate islands
 * 4) Assign p-value to each island using Stoffer-Liptak test
 * 5) Compute qvalues on islands p-values, filter by alpha
 *
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 * @param multiplier is used to
 * 1) Return broad peaks in case of broad modifications even for strict FDR settings
 * 2) Mitigate the problem when number of peaks for strict FDR is much bigger than for relaxed FDR
 */
fun SpanFitResults.getIslands(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    multiplier: Double = 1e2,
    cancellableState: CancellableState? = null
): List<Peak> {
    fitInfo.prepareScores()
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        // Check that we have information for requested chromosome
        if (chromosome.name in fitInfo.chromosomesSizes) {
            getChromosomeIslands(
                logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
                fitInfo.offsets(chromosome),
                chromosome,
                fdr,
                gap,
                multiplier
            )
        } else {
            SpanFitResults.LOG.debug(
                "NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}"
            )
            emptyList()
        }
    }
    return genomeQuery.get().flatMap { map[it] }
}


internal fun SpanFitResults.getChromosomeIslands(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    multiplier: Double,
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
    val resultIslands = filteredIslands.indices
        .filter { islandQValues[it] < fdr }
        .map { islandIndex ->
            val (i, j) = filteredIslands[islandIndex]
            val start = offsets[i]
            val end = if (j < offsets.size) offsets[j] else chromosome.length
            Peak(
                chromosome = chromosome,
                startOffset = start,
                endOffset = end,
                mlogpvalue = -log10(islandsPValues[islandIndex]),
                mlogqvalue = -log10(islandQValues[islandIndex]),
                // Value is either coverage of fold change
                value = fitInfo.score(ChromosomeRange(start, end, chromosome)),
                // Score should be proportional original q-value
                score = min(1000.0, -10 * log10(islandQValues[islandIndex])).toInt()
            )
        }
    Peak.LOG.debug(
        "$chromosome: islands result/filtered/candidate " +
                "${resultIslands.size}/${filteredIslands.size}/${candidateIslands.size}"
    )
    return resultIslands
}

/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * Using caches with weak values to avoid memory overflow.
 */
private val ISLANDS_STOFFER_LIPTAK_CACHE: Cache<Pair<SpanFitResults, Chromosome>, StofferLiptakTest> =
    CacheBuilder.newBuilder().weakValues().build()