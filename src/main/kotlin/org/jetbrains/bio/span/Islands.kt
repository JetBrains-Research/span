package org.jetbrains.bio.span

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.BenjaminiHochberg
import org.jetbrains.bio.statistics.hypothesis.StofferLiptakTest
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.floor
import kotlin.math.ln
import kotlin.math.log10

/**
 * The islands are called in five steps.
 *
 * 1) Estimate posterior probabilities
 * 2) Pick candidate bins with probability of H_0 <= fdr alpha
 * 3) Using gap merge bins into candidate islands
 * 4) Assign p-value to each island using Stoffer-Liptak test
 * 5) Compute qvalues on islands p-values
 *
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 */
fun SpanFitResults.getIslands(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    cancellableState: CancellableState? = null
): List<Peak> {
    val coverage = fitInfo.scoresDataFrame()
    if (coverage.isEmpty()) {
        SpanFitResults.LOG.debug("No coverage caches present, peak scores won't be computed.")
    }
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        // Check that we have information for requested chromosome
        if (chromosome.name in fitInfo.chromosomesSizes) {
            val f64LogNullMemberships = logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)
            val offsets = fitInfo.offsets(chromosome)
            getChromosomeIslands(
                f64LogNullMemberships, offsets, chromosome, fdr, gap, coverageDataFrame = coverage[chromosome]
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
    coverageDataFrame: DataFrame? = null
): List<Peak> {
    // Compute candidate bins and islands
    val candidateBins = BitterSet(logNullMemberships.size)
    val lnFdr = ln(fdr)
    0.until(candidateBins.size())
        .filter { logNullMemberships[it] <= lnFdr }
        .forEach(candidateBins::set)
    val candidateIslands = candidateBins.aggregate(gap)
    if (candidateIslands.isEmpty()) {
        return emptyList()
    }
    val nullMemberships = logNullMemberships.exp()
    val stofferLiptakTest = ISLANDS_STOFFER_LIPTAK_CACHE.get(Triple(this, chromosome, gap)) {
        StofferLiptakTest(nullMemberships.data)
    }

    val islandsPValues = ISLANDS_PVALUES_CACHE.get(Triple(this, chromosome, gap)) {
        // Apply Stoffer-Liptak test to correct dependence between consequent p-values
        F64Array(candidateIslands.size) { islandIndex ->
            val (i, j) = candidateIslands[islandIndex]
            stofferLiptakTest.combine(DoubleArray(j - i) { nullMemberships[it + i] })
        }
    }
    val islandQValues = ISLANDS_QVALUES_CACHE.get(Triple(this, chromosome, gap)) {
        BenjaminiHochberg.adjust(islandsPValues)
    }
    val minIslandQValues = islandQValues.min()
    val islands = candidateIslands.indices
        .filter { islandQValues[it] < fdr }
        .map { islandIndex ->
            val (i, j) = candidateIslands[islandIndex]
            val islandPValue = islandsPValues[islandIndex]
            val islandQValue = islandQValues[islandIndex]
            val islandStart = offsets[i]
            val islandEnd = if (j < offsets.size) offsets[j] else chromosome.length
            // Score should be proportional original q-value
            val score = floor(1000.0 * log10(islandQValue) / log10(minIslandQValues)).toInt()
            // Value is either coverage of fold change
            var value = 0.0
            if (coverageDataFrame != null) {
                if (coverageDataFrame.labels.size == 1 ||
                    coverageDataFrame.labels.all {
                        it.startsWith(SpanPeakCallingExperiment.TRACK_PREFIX)
                    }
                ) {
                    value = coverageDataFrame.partialMean(i, j)
                } else if (coverageDataFrame.labels.all {
                        it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK1_PREFIX) ||
                                it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK2_PREFIX)
                    }) {
                    val track1 = coverageDataFrame.partialMean(i, j, coverageDataFrame.labels
                        .filter {
                            it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK1_PREFIX)
                        })
                    val track2 = coverageDataFrame.partialMean(i, j, coverageDataFrame.labels
                        .filter {
                            it.startsWith(SpanDifferentialPeakCallingExperiment.TRACK2_PREFIX)
                        })
                    // Value if LogFC
                    value = if (track2 != 0.0) ln(track1) - ln(track2) else Double.MAX_VALUE
                } else {
                    Peak.LOG.debug("Failed to compute value for ${coverageDataFrame.labels}")
                }
            }
            Peak(
                chromosome, islandStart, islandEnd,
                mlogpvalue = -log10(islandPValue),
                mlogqvalue = -log10(islandQValue),
                value = value,
                score = score
            )
        }
    return islands
}


/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * Using caches with weak values to avoid memory overflow.
 */

private val ISLANDS_STOFFER_LIPTAK_CACHE: Cache<Triple<SpanFitResults, Chromosome, Int>, StofferLiptakTest> =
    CacheBuilder.newBuilder().weakValues().build()

private val ISLANDS_PVALUES_CACHE: Cache<Triple<SpanFitResults, Chromosome, Int>, F64Array> =
    CacheBuilder.newBuilder().weakValues().build()

private val ISLANDS_QVALUES_CACHE: Cache<Triple<SpanFitResults, Chromosome, Int>, F64Array> =
    CacheBuilder.newBuilder().weakValues().build()
