package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.statistics.hypothesis.StofferLiptakTest
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import java.lang.Integer.max
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.exp
import kotlin.math.ln
import kotlin.math.min


/**
 * The islands are called in five steps.
*
 * 1) Estimate posterior probabilities
 * 2) Pick candidate bins with relaxed posterior error probability, e.g. sqrt(fdr).
 * This mitigates the problem of wide marks peaks split on strong fdrs.
 * 3) Using gap merge bins into candidate islands.
 * 4) Assign p-value to each island using based on combined p-values for blocks of consequent enriched bins.
 *    Each block is assigned P as average posterior log error probability for bins in blocks.
 *    50% top significant blocks scores are aggregated using length-weighted average as P for island.
 * 5) Compute qvalues on islands p-values, filter by alpha.
*/
internal fun SpanFitResults.getIslands(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    summitsOnly: Boolean = false,
    cancellableState: CancellableState? = null
): List<Peak> {
    resetCounters()
    val progress = Progress { title = "Computing islands" }.bounded(genomeQuery.get().size.toLong())
    fitInfo.prepareScores()
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        val chromosomeIslands = getChromosomeIslands(chromosome, fdr, gap, summitsOnly)
        progress.report(1)
        chromosomeIslands
    }
    progress.done()
    Peak.LOG.debug(
        "Total candidate bins/candidate/result islands " +
                "${candidateBinsCounter.get()}/${candidateIslandsCounter.get()}/${resultIslandsCounter.get()}"
    )
    return genomeQuery.get().flatMap { map[it] }
}

internal fun SpanFitResults.getChromosomeIslands(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    summitsOnly: Boolean = false
): List<Peak> {
    // Check that we have information for requested chromosome
    val chromosomeIslands = if (chromosome.name in fitInfo.chromosomesSizes) {
        getChromosomeIslands(
            logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
            fitInfo.offsets(chromosome),
            chromosome,
            fdr,
            gap,
            summitsOnly
        )
    } else {
        SpanFitResults.LOG.debug(
            "NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}"
        )
        emptyList()
    }
    return chromosomeIslands
}

private val candidateBinsCounter = AtomicInteger()
private val candidateIslandsCounter = AtomicInteger()
private val resultIslandsCounter = AtomicInteger()

private fun resetCounters() {
    candidateBinsCounter.set(0)
    candidateIslandsCounter.set(0)
    resultIslandsCounter.set(0)
}

private fun SpanFitResults.getChromosomeIslands(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    summitsOnly: Boolean,
    relaxPower: Double = RELAX_POWER_DEFAULT,
): List<Peak> {
    /*
    * Compute candidate bins and islands with relaxed settings
    * Relaxed probability allows for:
    * 1) Return broad peaks in case of broad modifications even for strict FDR settings
    * 2) Mitigate the problem when number of peaks for strict FDR is much bigger than for relaxed FDR
    * */
    val logFdr = ln(fdr)
    val candidateBins = BitterSet(logNullMemberships.size).apply {
        0.until(size()).filter { logNullMemberships[it] <= min(ln(0.1), logFdr * relaxPower) }.forEach(::set)
    }
    val strictBins = BitterSet(logNullMemberships.size).apply {
        0.until(size()).filter { logNullMemberships[it] <= logFdr }.forEach(::set)
    }

    val candidateIslands = candidateBins.aggregate(gap).filter { (from, to) ->
        (from until to).any { logNullMemberships[it] <= logFdr }
    }
    if (candidateIslands.isEmpty()) {
        return emptyList()
    }

    // Init Stoffer-Liptak test to correct dependence between consequent p-values if required
    val slTest = if (summitsOnly) StofferLiptakTest(logNullMemberships.exp().toDoubleArray()) else null

    /**
     * We want two invariants from islands pvalues:
     * 1) The more strict FDR, the fewer peaks with smaller average length
     * 2) Peaks should not disappear when relaxing FDR
     *
     * Peak score is computed as length-weighted average p-value in its consequent enriched bins.
     */
    val islandsLogPeps = F64Array(candidateIslands.size) { islandIndex ->
        val blocks = strictBins.findConsequentBlocks(candidateIslands[islandIndex])
        val blocksLogPeps = blocks.map { (from, to) ->
            if (summitsOnly) {
                // Estimate only summits by sliding window in blocks - looking for summit with
                // the most significant Stoffer-Liptak test combined p-value
                val slidingWindowCombinedPs = SUMMIT_WINDOWS.flatMap { window ->
                    (max(0, from - window + 1)..min(to, offsets.size - window)).map { start ->
                        ln(slTest!!.combine(DoubleArray(window) { exp(logNullMemberships[start + it]) }))
                    }
                }
                // Minimum log p-value for any window size
                slidingWindowCombinedPs.minOrNull() ?: 1.0
            } else {
                // Average posterior log error probability for block
                KahanSum().apply {
                    (from until to).forEach { feed(logNullMemberships[it]) }
                }.result() / (to - from)
            }
        }

        /**
         * We want summary score to be robust wrt appending blocks of low significance,
         * so take into account top 50% blocks p-values, otherwise we'll get fdr-blinking peaks, i.e.
         * peaks which are present for stronger fdr, but missing for more relaxed settings
         * As MACS2 --broad use mean_from_value_length to take into account difference in blocks lengths
         */
        val sum = KahanSum()
        var l = 0
        blocks.zip(blocksLogPeps).sortedBy { it.second }.take(max(1, blocks.size / 2)).forEach { (b, p) ->
            sum += p * (b.toIndex - b.fromIndex)
            l += b.toIndex - b.fromIndex
        }
        sum.result() / l
    }

    // Filter result islands by Q values
    val islandsLogQValues = Fdr.qvalidate(islandsLogPeps, logResults = true)
    val resultIslandsIndexes = candidateIslands.indices.filter {
        islandsLogQValues[it] < logFdr
    }
    val resultIslands = resultIslandsIndexes.map { idx ->
        val (from, to) = candidateIslands[idx]
        val start = offsets[from]
        val end = if (to < offsets.size) offsets[to] else chromosome.length
        Peak(
            chromosome = chromosome,
            startOffset = start,
            endOffset = end,
            mlogpvalue = -islandsLogPeps[idx] / LOG_10,
            mlogqvalue = -islandsLogQValues[idx] / LOG_10,
            // Value is either coverage of fold change
            value = fitInfo.score(ChromosomeRange(start, end, chromosome)),
            // Score should be proportional original q-value
            score = min(1000.0, -10 * islandsLogQValues[idx] / LOG_10).toInt()
        )
    }
    Peak.LOG.debug(
        "$chromosome: candidate bins/candidate/result islands " +
                "${candidateBins.cardinality()}/${candidateIslands.size}/${resultIslands.size}"
    )
    candidateBinsCounter.addAndGet(candidateBins.cardinality())
    candidateIslandsCounter.addAndGet(candidateIslands.size)
    resultIslandsCounter.addAndGet(resultIslands.size)
    return resultIslands
}

const val RELAX_POWER_DEFAULT = 0.5
val SUMMIT_WINDOWS = intArrayOf(2, 3, 5)
