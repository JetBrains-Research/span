package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.dataframe.BitRange
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.span.coverage.CoverageScoresQuery.Companion.computeScales
import org.jetbrains.bio.span.fit.*
import org.jetbrains.bio.span.fit.SpanFitResults.Companion.LOG
import org.jetbrains.bio.span.fit.experimental.SpanRegrMixtureAnalyzeFitInformation
import org.jetbrains.bio.span.statistics.util.PoissonUtil
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.viktor.F64Array
import org.jetbrains.bio.viktor.KahanSum
import java.lang.Integer.max
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.ceil
import kotlin.math.ln
import kotlin.math.min


object SpanModelToPeaks {

    /**
     * Main method to compute peaks from model.
     *
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
    fun computeChromosomePeaks(
        spanFitResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        fdr: Double,
        gap: Int,
        cancellableState: CancellableState? = null
    ): List<Peak> {
        resetCounters()
        val progress = Progress { title = "Computing peaks fdr=$fdr gap=$gap" }.bounded(genomeQuery.get().size.toLong())
        spanFitResults.fitInfo.prepareScores()
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            cancellableState?.checkCanceled()
            val chromosomePeaks = computeChromosomePeaks(spanFitResults, chromosome, fdr, gap)
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


    fun computeChromosomePeaks(
        spanFitResults: SpanFitResults,
        chromosome: Chromosome,
        fdr: Double,
        gap: Int,
    ): List<Peak> {
        // Check that we have information for requested chromosome
        val chromosomeIslands = if (chromosome.name in spanFitResults.fitInfo.chromosomesSizes) {
            getChromosomePeaks(
                chromosome,
                spanFitResults.fitInfo,
                spanFitResults.logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL),
                spanFitResults.fitInfo.offsets(chromosome),
                fdr,
                gap,
            )
        } else {
            LOG.debug("NO peaks information for chromosome: ${chromosome.name} in fitInfo ${spanFitResults.fitInfo.build}")
            emptyList()
        }
        return chromosomeIslands
    }

    private fun getChromosomePeaks(
        chromosome: Chromosome,
        fitInfo: SpanFitInformation,
        logNullMemberships: F64Array,
        offsets: IntArray,
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
        val treatmentCoverage: Coverage?
        val controlCoverage: Coverage?
        var treatmentScale = 1.0
        var controlScale = 1.0
        when (fitInfo) {
            is SpanAnalyzeFitInformation -> {
                val treatmentReads = fitInfo.scoreQueries!!.single().treatmentReads
                treatmentCoverage = if (treatmentReads.isAccessible()) treatmentReads.get() else null
                val controlReads = fitInfo.scoreQueries!!.single().controlReads
                controlCoverage = if (controlReads?.isAccessible() == true) controlReads.get() else null
            }
            is SpanCompareFitInformation -> {
                val treatmentReads = fitInfo.scoreQueries1!!.single().treatmentReads
                treatmentCoverage = if (treatmentReads.isAccessible()) treatmentReads.get() else null
                val controlReads = fitInfo.scoreQueries2!!.single().treatmentReads
                controlCoverage = if (controlReads.isAccessible()) controlReads.get() else null
            }
            is SpanRegrMixtureAnalyzeFitInformation -> {
                val treatmentReads = fitInfo.scoreQuery!!.treatmentReads
                treatmentCoverage = if (treatmentReads.isAccessible()) treatmentReads.get() else null
                controlCoverage = null
            }
            else -> {
                throw IllegalStateException("Incorrect fitInfo: ${fitInfo.javaClass.name}")
            }
        }
        if (treatmentCoverage != null && controlCoverage != null) {
            val scales = computeScales(fitInfo.genomeQuery(), treatmentCoverage, controlCoverage)!!
            treatmentScale = scales.first
            controlScale = scales.second
        }
        val islandsLogPs = F64Array(candidateIslands.size) { islandIndex ->
            val blocks = strictBins.findConsequentBlocks(candidateIslands[islandIndex])
            val blocksLogPs = blocks.map { (from, to) ->
                if (treatmentCoverage != null && controlCoverage != null) {
                    // Estimate enrichment vs local coverage in control track
                    val start = offsets[from]
                    val end = if (to < offsets.size) offsets[to] else chromosome.length
                    val chromosomeRange = ChromosomeRange(start, end, chromosome)
                    // Scaling down by (to - from) allows to align p-values,
                    // but results in less peaks in low frip conditions
                    val lengthScale = 10.0 / (to - from)
                    val peakTreatment =
                        treatmentCoverage.getBothStrandsCoverage(chromosomeRange) * treatmentScale / lengthScale
                    val peakControl =
                        controlCoverage.getBothStrandsCoverage(chromosomeRange) * controlScale / lengthScale
                    PoissonUtil.logPoissonCdf(ceil(peakTreatment).toInt() + pseudoCount, peakControl + pseudoCount)
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

    /**
     * We want summary score to be robust wrt appending blocks of low significance,
     * so take into account top 50% blocks p-values, otherwise we'll get fdr-blinking peaks, i.e.
     * peaks which are present for stronger fdr, but missing for more relaxed settings
     * As MACS2 --broad use mean_from_value_length to take into account difference in blocks lengths
     */
    private fun islandsLengthWeightedScores(
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

    private fun relaxedLogFdr(logFdr: Double, relaxPower: Double) = min(ln(0.1), logFdr * relaxPower)

    private const val RELAX_POWER_DEFAULT = 0.5

    private val LOG_10 = ln(10.0)

    private val candidateBinsCounter = AtomicInteger()
    private val candidateIslandsCounter = AtomicInteger()
    private val resultIslandsCounter = AtomicInteger()

    private fun resetCounters() {
        candidateBinsCounter.set(0)
        candidateIslandsCounter.set(0)
        resultIslandsCounter.set(0)
    }

}