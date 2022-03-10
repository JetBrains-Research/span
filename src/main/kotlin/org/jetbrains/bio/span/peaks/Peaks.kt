package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.dataframe.BitRange
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.span.coverage.CoverageScoresQuery.Companion.computeScales
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanCompareFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.fit.SpanFitResults.Companion.LOG
import org.jetbrains.bio.span.fit.SpanModelFitExperiment
import org.jetbrains.bio.span.fit.experimental.SpanRegrMixtureAnalyzeFitInformation
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import org.jetbrains.bio.viktor.logAddExp
import java.lang.Integer.max
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.ln
import kotlin.math.min


/**
 * 1) Estimate posterior probabilities
 * 2) Pick candidate bins with relaxed posterior error probability, e.g. sqrt(fdr).
 * This mitigates the problem of wide marks peaks split on strong fdrs.
 * 3) Using gap merge bins into candidate islands.
 * 4) Assign p-value to each island using based on combined p-values for blocks of consequent enriched bins.
 *    In case when control track is present, we use Poisson CDF to estimate log P value,
 *    otherwise an average log PEP (posterior error probability) for bins in blocks is used.
 * 5) 50% top significant blocks scores are aggregated using length-weighted average as P for island.
 * 6) Compute qvalues by islands p-values, filter by alpha.
 */
fun SpanFitResults.getPeaks(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    cancellableState: CancellableState? = null
): List<Peak> {
    resetCounters()
    val progress = Progress { title = "Computing peaks fdr=$fdr gap=$gap" }.bounded(genomeQuery.get().size.toLong())
    fitInfo.prepareScores()
    val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
        cancellableState?.checkCanceled()
        val chromosomePeaks = getPeaks(chromosome, fdr, gap)
        progress.report(1)
        chromosomePeaks
    }
    progress.done()
    Peak.LOG.debug(
        "Total candidate bins/candidate/result islands " +
                "${candidateBinsCounter.get()}/${candidateIslandsCounter.get()}/${resultIslandsCounter.get()}"
    )
    return genomeQuery.get().flatMap { map[it] }
}

fun SpanFitResults.getPeaks(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
): List<Peak> {
    // Check that we have information for requested chromosome
    val chromosomeIslands = if (chromosome.name in fitInfo.chromosomesSizes) {
        getChromosomePeaks(
            logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
            fitInfo.offsets(chromosome),
            chromosome,
            fdr,
            gap,
        )
    } else {
        LOG.debug("NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}")
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

private fun SpanFitResults.getChromosomePeaks(
    logNullMemberships: F64Array,
    offsets: IntArray,
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    relaxPower: Double = RELAX_POWER_DEFAULT,
    pseudoCount: Int = 1
): List<Peak> {
    // Compute candidate bins and islands with relaxed settings
    // Relaxed probability allows for:
    // 1) Return broad peaks in case of broad modifications even for strict FDR settings
    // 2) Mitigate the problem when number of peaks for strict FDR is much bigger than for relaxed FDR
    val logFdr = ln(fdr)
    val strictBins = BitterSet(logNullMemberships.size) { logNullMemberships[it] <= logFdr }
    val relaxedLogFdr = relaxedLogFdr(logFdr, relaxPower)
    val candidateBins = BitterSet(logNullMemberships.size).apply {
        0.until(size()).filter { logNullMemberships[it] <= relaxedLogFdr }.forEach(::set)
    }
    val candidateIslands = candidateBins.aggregate(gap).filter { (from, to) ->
        (from until to).any { logNullMemberships[it] <= logFdr }
    }
    if (candidateIslands.isEmpty()) {
        return emptyList()
    }
    // We want two invariants from islands pvalues:
    // 1) The more strict FDR, the fewer peaks with smaller average length
    // 2) Peaks should not disappear when relaxing FDR
    // Peak score is computed as length-weighted average p-value in its consequent enriched bins.
    val treatmentCoverage: Coverage
    val controlCoverage: Coverage?
    var treatmentScale = 1.0
    var controlScale = 1.0
    if (fitInfo is SpanAnalyzeFitInformation) {
        treatmentCoverage = fitInfo.scoreQueries!!.single().treatmentReads.get()
        controlCoverage = fitInfo.scoreQueries!!.single().controlReads?.get()
        if (controlCoverage != null) {
            val scales = computeScales(fitInfo.genomeQuery(), treatmentCoverage, controlCoverage)!!
            treatmentScale = scales.first
            controlScale = scales.second
        }
    } else if (fitInfo is SpanCompareFitInformation) {
        treatmentCoverage = fitInfo.scoreQueries1!!.single().treatmentReads.get()
        controlCoverage = fitInfo.scoreQueries2!!.single().treatmentReads.get()
        val scales = computeScales(fitInfo.genomeQuery(), treatmentCoverage, controlCoverage)!!
        treatmentScale = scales.first
        controlScale = scales.second
    } else if (fitInfo is SpanRegrMixtureAnalyzeFitInformation) {
        treatmentCoverage = fitInfo.scoreQuery!!.treatmentReads.get()
        controlCoverage = null
    } else {
        throw IllegalStateException("Incorrect fitInfo: ${fitInfo.javaClass.name}")
    }
    val islandsLogPs = F64Array(candidateIslands.size) { islandIndex ->
        val blocks = strictBins.findConsequentBlocks(candidateIslands[islandIndex])
        val blocksLogPs = blocks.map { (from, to) ->
            if (controlCoverage != null) {
                // Estimate enrichment vs local coverage in control track
                val start = offsets[from]
                val end = if (to < offsets.size) offsets[to] else chromosome.length
                val chromosomeRange = ChromosomeRange(start, end, chromosome)
                // Scaling down by (to - from) allows to align p-values,
                // but results in less peaks in low frip conditions
                val lengthScale = 10.0 / (to - from)
                val peakTreatment =
                    treatmentCoverage.getBothStrandsCoverage(chromosomeRange) * treatmentScale / lengthScale
                val peakControl = controlCoverage.getBothStrandsCoverage(chromosomeRange) * controlScale / lengthScale
                log10PoissonCdfQLargeLambda(ceil(peakTreatment).toInt() + pseudoCount, peakControl + pseudoCount)
            } else {
                // Average posterior log error probability for block
                KahanSum().apply {
                    (from until to).forEach { feed(logNullMemberships[it]) }
                }.result() / (to - from)
            }
        }
        islandsLengthWeightedScores(blocks, blocksLogPs)
    }

    // Filter result islands by Q values
    val islandsLogQValues = Fdr.qvalidate(islandsLogPs, logResults = true)
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
            mlogpvalue = -islandsLogPs[idx] / LOG_10,
            mlogqvalue = -islandsLogQValues[idx] / LOG_10,
            // Value is either coverage of fold change
            value = fitInfo.score(ChromosomeRange(start, end, chromosome)),
            // Score should be proportional original q-value
            score = min(1000.0, -10 * islandsLogQValues[idx] / LOG_10).toInt()
        )
    }

    Peak.LOG.debug(
        "$chromosome: candidate bins / candidate/result islands " +
                "${candidateBins.cardinality()} / ${candidateIslands.size}/${resultIslands.size}"
    )
    candidateBinsCounter.addAndGet(candidateBins.cardinality())
    candidateIslandsCounter.addAndGet(candidateIslands.size)
    resultIslandsCounter.addAndGet(resultIslands.size)
    return resultIslands
}

fun relaxedLogFdr(logFdr: Double, relaxPower: Double) = min(ln(0.1), logFdr * relaxPower)

const val RELAX_POWER_DEFAULT = 0.5

val LOG_10 = ln(10.0)

/**
 * We want summary score to be robust wrt appending blocks of low significance,
 * so take into account top 50% blocks p-values, otherwise we'll get fdr-blinking peaks, i.e.
 * peaks which are present for stronger fdr, but missing for more relaxed settings
 * As MACS2 --broad use mean_from_value_length to take into account difference in blocks lengths
 */
fun islandsLengthWeightedScores(
    blocks: List<BitRange>,
    scores: List<Double>,
    fraction: Double = 0.5
): Double {
    require(blocks.size == scores.size) { "Different lengths of blocks and scores lists" }
    val sum = KahanSum()
    var l = 0
    blocks.zip(scores).sortedBy { it.second }.take(max(1, (blocks.size * fraction).toInt())).forEach { (b, p) ->
        sum += p * (b.toIndex - b.fromIndex)
        l += b.toIndex - b.fromIndex
    }
    return sum.result() / l
}

/**
 * Slower Poisson CDF evaluater for upper tail which allow calculation in log space.
 * Better for the pvalue < 10^-310.
 * @param k observation
 * @param lbd: Lambda
 * @return log(pvalue)
 *
 * See MACS2 sources Prob.pyx for original source code:
 * ret = -lambda + \ln( \sum_{i=k+1}^{\inf} {lambda^i/i!} = -lambda + \ln( sum{ exp{ln(F)} } ), where F=lambda^m/m!
 * \ln{F(m)} = m*ln{lambda} - \sum_{x=1}^{m}\ln(x)
 * Calculate \ln( sum{exp{N} ) by logspace_add function
 */
fun log10PoissonCdfQLargeLambda(k: Int, lbd: Double, maxM: Int = 10000, epsilon: Double = 1e-5): Double {
    require(lbd > 0) {
        "Lambda should be > 0, got $lbd"
    }
    var residue: Double
    var logX: Double
    val lnLbd = ln(lbd)
    // first residue
    val m = k + 1
    var sumLns = 0.0 // TODO[shpynov] this may be tabulated
    for (i in 1 until m + 1) {
        sumLns += ln(i.toDouble())
    }
    logX = m * lnLbd - sumLns
    residue = logX
    var logy: Double
    for (i in m + 1..maxM) { // Limit
        logy = logX + lnLbd - ln(i.toDouble())
        val preResidue = residue
        residue = preResidue logAddExp logy
        if (abs(preResidue - residue) < epsilon)
            break
        logX = logy
    }
    return residue - lbd
}