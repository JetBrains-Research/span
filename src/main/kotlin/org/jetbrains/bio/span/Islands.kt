package org.jetbrains.bio.span

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.viktor.F64Array
import kotlin.math.ln
import kotlin.math.log10
import kotlin.math.min

/**
 * The islands are called in five steps.
 *
 * 1) Estimate posterior probabilities
 * 2) Pick bins with probability of H_0 <= nullProbabilityThreshold
 * 3) Using gap merge bins into islands (before q-values)
 * 4) Assign score to each island as median(log(ps)) * log(length)
 * 5) Compute qvalues on scores
 *
 * @param offsets All the data is binarized, so offsets stores positions in base pair.
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 * @param nullProbabilityThreshold threshold to collect candidate islands
 * @param coverageDataFrame is used to compute either coverage or log fold change.
 */
internal fun SpanFitResults.getChromosomeIslands(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    nullProbabilityThreshold: Double,
    coverageDataFrame: DataFrame? = null
): List<Peak> {
    val logNullProbabilityThreshold = ln(nullProbabilityThreshold)
    val candidateBins = BitterSet(logNullMemberships.size)
    0.until(candidateBins.size())
        .filter { logNullMemberships[it] <= logNullProbabilityThreshold }
        .forEach(candidateBins::set)
    val candidateIslands = candidateBins.aggregate(gap)
    if (candidateIslands.isEmpty()) {
        return emptyList()
    }
    // Assign p-value like scores to merged peaks using SICER inspired scheme
    val islandsLogNullMemberships = F64Array(candidateIslands.size) { islandIndex ->
        // SICER scheme sum(-log(ps)) leads to huge bias towards long islands
        // Use median * log(length)
        val (i, j) = candidateIslands[islandIndex]
        val medianLogNullMembership = (i until j)
            .map { logNullMemberships[it] }
            .filter { it <= logNullProbabilityThreshold }
            .toDoubleArray().median()
        medianLogNullMembership * ln((j - i).toDouble())
    }
    val islandQValues = islandsQValuesCache.get(Triple(this, chromosome, gap)) {
        Fdr.qvalidate(islandsLogNullMemberships)
    }
    val islands = candidateIslands.indices
        .filter { islandQValues[it] < fdr }
        .map { islandIndex ->
            val (i, j) = candidateIslands[islandIndex]
            val islandLogNullMembership = islandsLogNullMemberships[islandIndex]
            val islandQValue = islandQValues[islandIndex]
            val islandStart = offsets[i]
            val islandEnd = if (j < offsets.size) offsets[j] else chromosome.length
            // Score should be proportional to length of peak and average original q-value
            val score = min(1000.0, -log10(islandQValue))
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
                mlogpvalue = -islandLogNullMembership,
                mlogqvalue = -log10(islandQValue),
                value = value,
                score = score.toInt()
            )
        }
    return islands
}

fun SpanFitResults.getIslands(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    nullProbabilityThreshold: Double,
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
            val f64LogNullMemberships =
                logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)
            val offsets = fitInfo.offsets(chromosome)
            getChromosomeIslands(f64LogNullMemberships, offsets, chromosome,
                fdr, gap, nullProbabilityThreshold, coverage[chromosome])
        } else {
            SpanFitResults.LOG.debug("NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}")
            emptyList()
        }
    }
    return genomeQuery.get().flatMap { map[it] }
}

const val SPAN_ISLANDS_DEFAULT_NULL_PROBABILITY = 0.2


/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * So Q-values estimation is superfluous for each parameters combination.
 * Use cache with weak values to avoid memory overflow.
 */
private val islandsQValuesCache: Cache<Triple<SpanFitResults, Chromosome, Int>, F64Array> =
    CacheBuilder.newBuilder().weakValues().build()
