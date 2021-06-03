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
 * The islands are called in three steps.
 *
 * Estimate posterior probabilities
 * Pick enriched bins with q-value <= FDR
 * Extend enriched bins into islands with given GAP and other bins, with logMembership <= FDR
 * Assign score to each island proportional to length and peak fdr
 *
 * @param offsets All the data is binarized, so offsets stores positions in base pair.
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 * @param coverageDataFrame is used to compute either coverage or log fold change.
 */
fun SpanFitResults.getChromosomeIslands(
    logNullMemberships: F64Array,
    qvalues: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    coverageDataFrame: DataFrame? = null
): List<Peak> {
    // Check that we have information for requested chromosome
    if (chromosome.name in fitInfo.chromosomesSizes) {
        // Mark enriched bins by pvalue
        val enrichedBins = BitterSet(logNullMemberships.size)
        0.until(enrichedBins.size()).filter { logNullMemberships[it] < ln(fdr) }.forEach(enrichedBins::set)

        // Aggregate enriched islands, and pick those, which contain enriched bin by fdr
        return enrichedBins.aggregate(gap)
            .filter { (i, j) ->
                (i until j).any { qvalues[it] < fdr }
            }.map { (i, j) ->
                val passedFDR = (i until j).filter { qvalues[it] < fdr }
                val pvalueLogMedian = DoubleArray(passedFDR.size) { logNullMemberships[passedFDR[it]] }.median()
                val qvalueMedian = DoubleArray(passedFDR.size) { qvalues[passedFDR[it]] }.median()
                val start = offsets[i]
                val end = if (j < offsets.size) offsets[j] else chromosome.length
                // Score should be proportional to length of peak and median of original q-value
                val score = min(1000.0, (-log10(qvalueMedian) * (1 + ln((end - start).toDouble()))))
                // Value is either coverage of fold change
                var value = 0.0
                if (coverageDataFrame != null) {
                    if (coverageDataFrame.labels.size == 1 ||
                        coverageDataFrame.labels.all { it.startsWith(SpanPeakCallingExperiment.TRACK_PREFIX) }
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
                    chromosome, start, end,
                    mlogpvalue = -pvalueLogMedian,
                    mlogqvalue = -log10(qvalueMedian),
                    value = value,
                    score = score.toInt()
                )
            }
    } else {
        SpanFitResults.LOG.debug("NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}")
        return emptyList()
    }
}

/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * So Q-values estimation is superfluous for each parameters combination.
 * Use cache with weak values to avoid memory overflow.
 */
private val f64QValuesCache: Cache<Pair<SpanFitResults, Chromosome>, F64Array> =
    CacheBuilder.newBuilder().weakValues().build()

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
            val f64LogNullMemberships =
                logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)
            val f64QValues = f64QValuesCache.get(this to chromosome) {
                Fdr.qvalidate(f64LogNullMemberships)
            }
            val offsets = fitInfo.offsets(chromosome)
            getChromosomeIslands(f64LogNullMemberships, f64QValues, offsets, chromosome, fdr, gap, coverage[chromosome])
        } else {
            SpanFitResults.LOG.debug("NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}")
            emptyList()
        }
    }
    return genomeQuery.get().flatMap { map[it] }
}
