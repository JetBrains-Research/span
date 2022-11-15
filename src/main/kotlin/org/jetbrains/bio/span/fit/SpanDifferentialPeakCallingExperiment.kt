package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.span.coverage.BinnedNormalizedCoverageQuery
import org.jetbrains.bio.span.peaks.ModelToPeaks
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.statistics.hmm.ConstrainedNBZHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.MultiLabels
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * Corresponds to Span `compare` invocation.
 *
 * The treatment-control pairs are split into two sets that are to compare.
 *
 * For each treatment-control pair, we compute binned normalized coverage [BinnedNormalizedCoverageQuery].
 * These coverages are used as the input for a five-state multidimensional negative binomial HMM.
 * For each dimension `d`, there are two negative binomial distributions, low_d and high_d.
 * - ZERO state corresponds to zero emissions for all dimensions
 * - LOW state employs `low_d` emission for each dimension `d`
 * - HIGH state employs `high_d` emission for each dimension `d`
 * - INCREASED state employs `low_d` emission for each dimension `d` from the first set and `high_d` for the second set
 * - DECREASED state employs `high_d` emission for each dimension `d` from the first set and `low_d` for the second set
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment private constructor(
    fitInformation: SpanCompareFitInformation,
    threshold: Double,
    maxIterations: Int
) : SpanModelFitExperiment<ConstrainedNBZHMM, SpanCompareFitInformation, ZLHID>(
    fitInformation,
    ConstrainedNBZHMM.fitter(fitInformation.data1.size, fitInformation.data2.size),
    ConstrainedNBZHMM::class.java,
    ZLHID.values(), NullHypothesis.of(ZLHID.same()),
    threshold = threshold,
    maxIterations = maxIterations
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

    fun computeDirectedDifferencePeaks(
        fdr: Double,
        gap: Int
    ): Pair<List<Peak>, List<Peak>> {
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            ModelToPeaks.computeChromosomePeaks(results, chromosome, fdr, gap, true)
        }
        val highLow = arrayListOf<Peak>()
        val lowHigh = arrayListOf<Peak>()
        genomeQuery.get().forEach { chromosome ->
            val states = getStatesDataFrame(chromosome)
            map[chromosome].forEach {
                if (states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.D.name) >
                    states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.I.name)
                ) {
                    highLow.add(it)
                } else {
                    lowHigh.add(it)
                }
            }
        }
        return highLow to lowHigh
    }


    companion object {
        private const val TRACK1_PREFIX = "track1"
        private const val TRACK2_PREFIX = "track2"

        /**
         * Creates experiment for model-based comparison of binned coverage tracks for given queries.
         *
         * @return experiment [SpanDifferentialPeakCallingExperiment]
         */
        fun getExperiment(
            genomeQuery: GenomeQuery,
            paths1: List<SpanDataPaths>,
            paths2: List<SpanDataPaths>,
            bin: Int,
            fragment: Fragment,
            unique: Boolean,
            threshold: Double,
            maxIterations: Int
        ): SpanDifferentialPeakCallingExperiment {
            require(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            val fitInformation = SpanCompareFitInformation.effective(
                genomeQuery,
                paths1, paths2,
                MultiLabels.generate(TRACK1_PREFIX, paths1.size).toList(),
                MultiLabels.generate(TRACK2_PREFIX, paths2.size).toList(),
                fragment, unique, bin
            )
            return SpanDifferentialPeakCallingExperiment(fitInformation, threshold, maxIterations)
        }
    }
}