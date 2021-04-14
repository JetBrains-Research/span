package org.jetbrains.bio.span

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.dataframe.BitterSet
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.*
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.statistics.hypothesis.Fdr
import org.jetbrains.bio.util.CancellableState
import org.jetbrains.bio.util.bufferedWriter
import org.jetbrains.bio.viktor.F64Array
import java.nio.file.Path
import kotlin.math.ln
import kotlin.math.log10
import kotlin.math.min

/**
 * The islands are called in three steps.
 *
 * Estimate posterior probabilities
 * Pick bins with probability of signal >= 0.8 (THIS is very important for sensitivity / peak lengths)
 * Using gap merge bins into islands (before q-values)
 * Assign score to each island as median(log(ps)) * log(length)
 * (SICER scheme as sum(-log(ps)) leads to huge bias towards long islands
 *  and inefficient qvalue filtering
 *  THIS can influence further qvalues estimation)
 * Compute qvalues on scores
 *
 * @param offsets All the data is binarized, so offsets stores positions in base pair.
 * @param fdr is used to limit False Discovery Rate at given level.
 * @param gap enriched bins yielded after FDR control are merged if distance is less or equal than gap.
 * @param coverageDataFrame is used to compute either coverage or log fold change.
 */
fun SpanFitResults.getChromosomeIslands(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    coverageDataFrame: DataFrame? = null,
    nullProbabilityThreshold: Double = 0.2
): List<Peak> {
    // Check that we have information for requested chromosome
    if (chromosome.name in fitInfo.chromosomesSizes) {
        val f64LogNullMemberships = logNullMemberships[chromosome.name]!!.f64Array(SpanModelFitExperiment.NULL)
        val candidateBins = BitterSet(f64LogNullMemberships.size)
        0.until(candidateBins.size())
            .filter { f64LogNullMemberships[it] <= ln(nullProbabilityThreshold) }
            .forEach(candidateBins::set)
        val candidateIslands = candidateBins.aggregate(gap)
        if (candidateIslands.isEmpty()) {
            return emptyList()
        }
        val islandsLogNullMemberships = F64Array(candidateIslands.size) { islandIndex ->
            val (i, j) = candidateIslands[islandIndex]
            val bins = (i until j).filter { f64LogNullMemberships[it] <= nullProbabilityThreshold }
            DoubleArray(bins.size) { f64LogNullMemberships[bins[it]] }.median() * ln((j - i).toDouble())
            // This score scheme creates huge bias towards long islands with moderate enrichment level,
            // which leads to quite small reduce in peaks number for stringent FDR 0.1 = 50K to 1E-20 = 30K
            // DoubleArray(bins.size) { f64LogNullMemberships[bins[it]] }.sum()
        }
        val islandQValues = islandsQValuesCache.get(Triple(this, chromosome, gap)) {
            Fdr.qvalidate(islandsLogNullMemberships)
        }
        val offsets = fitInfo.offsets(chromosome)
        return candidateIslands.indices
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
                        }) {
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
                Peak(chromosome, islandStart, islandEnd,
                    mlogpvalue = -islandLogNullMembership,
                    mlogqvalue = -log10(islandQValue),
                    value = value,
                    score = score.toInt())
            }
    } else {
        SpanFitResults.LOG.debug("NO peaks information for chromosome: ${chromosome.name} in fitInfo ${fitInfo.build}")
        return emptyList()
    }
}

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
        getChromosomeIslands(chromosome, fdr, gap, coverage[chromosome])
    }
    return genomeQuery.get().flatMap { map[it] }
}

fun saveIslands(peaks: List<Peak>, path: Path, peakName: String = "") {
    Peak.LOG.debug("""FORMAT $path:
chromosome, start, end, name, score, strand, coverage/foldchange, -log(pvalue), -log(qvalue)""")
    CSVFormat.TDF.print(path.bufferedWriter()).use { printer ->
        peaks.sorted().forEachIndexed { i, peak ->
            /* See MACS2 output format for details https://github.com/taoliu/MACS/ */
            printer.printRecord(
                peak.chromosome.name,
                peak.range.startOffset.toString(),
                peak.range.endOffset.toString(),
                "${if (peakName.isNotEmpty()) peakName else "peak"}_${i + 1}",
                peak.score.toString(),
                ".",
                peak.value.toString(),
                peak.mlogpvalue.toString(),
                peak.mlogqvalue.toString())
        }
    }
}

/**
 * During SPAN models optimizations we iterate over different FDR and GAPs parameters,
 * So Q-values estimation is superfluous for each parameters combination.
 * Use cache with weak values to avoid memory overflow.
 */
private val islandsQValuesCache: Cache<Triple<SpanFitResults, Chromosome, Int>, F64Array> =
    CacheBuilder.newBuilder().weakValues().build()
