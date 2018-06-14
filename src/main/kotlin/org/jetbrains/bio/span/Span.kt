package org.jetbrains.bio.span

import com.google.common.annotations.VisibleForTesting
import org.apache.log4j.Logger
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.experiments.fit.ModelFitExperiment
import org.jetbrains.bio.experiments.fit.SpanDifferentialPeakCallingExperiment
import org.jetbrains.bio.experiments.fit.SpanPeakCallingExperiment
import org.jetbrains.bio.experiments.tuning.SPAN
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.InputQuery
import org.jetbrains.bio.statistics.*
import org.jetbrains.bio.statistics.data.DataFrame
import org.jetbrains.bio.statistics.hmm.MLAbstractHMM
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.state.ZHL
import org.jetbrains.bio.statistics.state.ZLHID


/**
 * [Span] (Semi-supervised Peak Analyzer) is a tool for analyzing and comparing ChIP-Seq data.
 * Both procedures rely on the Zero Inflated Negative Binomial Restricted Algorithm.
 *
 * It is implemented as [ModelFitExperiment] with different [ClassificationModel] models.
 *
 * Enrichment
 * - States: [ZHL]
 * - Single replicate: [MLFreeNBHMM] zero-inflated HMM with univariate Negative Binomial emissions
 * - Multi replicates: [MLConstrainedNBHMM] zero-inflated HMM with multidimensional Negative Binomial emissions
 *
 * Difference
 * - States: [ZLHID]
 * - Any number of replicates: [MLConstrainedNBHMM]
 */
object Span {
    private val LOG = Logger.getLogger(Span::class.java)
    const val BIN = SPAN.DEFAULT_BIN
    const val FDR = SPAN.DEFAULT_FDR
    const val GAP = SPAN.DEFAULT_GAP


    /**
     * Creates experiment for model-based enrichment of binned coverage tracks (e.g. ChIP-seq tracks)
     * for given number of [queries].
     * Not restricted for single query and constrained for multiple queries.
     *
     * @return experiment [SpanPeakCallingExperiment]
     */
    fun getPeakCallingExperiment(genomeQuery: GenomeQuery,
                                 queries: List<InputQuery<Coverage>>,
                                 bin: Int):
            SpanPeakCallingExperiment<out MLAbstractHMM, ZHL> {
        check(queries.isNotEmpty()) { "No data" }
        return if (queries.size == 1) {
            SpanPeakCallingExperiment(genomeQuery, queries.first(),
                    semanticCheck(MLFreeNBHMM.fitter()),
                    MLFreeNBHMM::class.java,
                    bin, ZHL.values(), NullHypothesis.of(ZHL.Z, ZHL.L))
        } else {
            SpanPeakCallingExperiment(genomeQuery, queries,
                    semanticCheck(MLConstrainedNBHMM.fitter(queries.size), queries.size),
                    MLConstrainedNBHMM::class.java,
                    bin, ZHL.values(), NullHypothesis.of(ZHL.Z, ZHL.L))
        }
    }

    /**
     * Creates experiment for model-based comparison of binned coverage tracks for given queries.
     * Not restricted for single query and constrained for multiple queries.
     *
     * @return experiment [SpanDifferentialPeakCallingExperiment]
     */
    fun getDifferentialPeakCallingExperiment(genomeQuery: GenomeQuery,
                                             queries1: List<InputQuery<Coverage>>,
                                             queries2: List<InputQuery<Coverage>>,
                                             bin: Int):
            SpanDifferentialPeakCallingExperiment<*, ZLHID> {
        check(queries1.isNotEmpty() && queries2.isNotEmpty()) { "No data" }
        return if (queries1.size == 1 && queries2.size == 1) {
            SpanDifferentialPeakCallingExperiment(genomeQuery, queries1.first(), queries2.first(), bin,
                    semanticCheck(MLConstrainedNBHMM.fitter(1, 1), 1, 1),
                    MLConstrainedNBHMM::class.java,
                    ZLHID.values(), NullHypothesis.of(ZLHID.same()))
        } else {
            SpanDifferentialPeakCallingExperiment(genomeQuery, queries1, queries2, bin,
                    semanticCheck(MLConstrainedNBHMM.fitter(queries1.size, queries2.size),
                            queries1.size, queries2.size),
                    MLConstrainedNBHMM::class.java,
                    ZLHID.values(), NullHypothesis.of(ZLHID.same()))
        }
    }

    private fun semanticCheck(fitter: Fitter<MLConstrainedNBHMM>, tracks1: Int, tracks2: Int): Fitter<MLConstrainedNBHMM> {
        return object : Fitter<MLConstrainedNBHMM> by fitter {

            override fun fit(preprocessed: Preprocessed<DataFrame>,
                             title: String,
                             threshold: Double,
                             maxIter: Int): MLConstrainedNBHMM {
                val model = fitter.fit(preprocessed, title, threshold, maxIter)
                stateFlip(model, tracks1, tracks2)
                return model
            }

            override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                             title: String,
                             threshold: Double,
                             maxIter: Int): MLConstrainedNBHMM {
                val model = fitter.fit(preprocessed, title, threshold, maxIter)
                stateFlip(model, tracks1, tracks2)
                return model
            }
        }
    }

    /**
     * Flip states in case when states with HIGH get lower mean than LOW
     */
    @VisibleForTesting
    internal fun stateFlip(model: MLConstrainedNBHMM, tracks1: Int, tracks2: Int) {
        val means = model.means
        val ps = model.successProbabilities
        val switchNeeded1 = (0 until tracks1).filter { means[it] > means[it + tracks1] && ps[it] > ps[it + tracks1] }
        val switchNotNeeded1 = (0 until tracks1).filter { means[it] < means[it + tracks1] && ps[it] < ps[it + tracks1] }
        val switchNeeded2 = (2 * tracks1 until 2 * tracks1 + tracks2)
                .filter { means[it] > means[it + tracks2] && ps[it] > ps[it + tracks2] }
        val switchNotNeeded2 = (2 * tracks1 until 2 * tracks1 + tracks2)
                .filter { means[it] < means[it + tracks2] && ps[it] < ps[it + tracks2] }
        if (switchNeeded1.isNotEmpty() && switchNotNeeded1.isNotEmpty()) {
            LOG.error("Irrecoverable fitting error")
            LOG.error("means: " + means.toString())
            LOG.error("ps: " + ps.toString())
            LOG.error("track(s) " + switchNeeded1.joinToString(transform = Int::toString)
                    + " contradict track(s) " + switchNotNeeded1.joinToString(transform = Int::toString))
            throw IllegalStateException("Irrecoverable fitting error")
        }
        if (switchNeeded2.isNotEmpty() && switchNotNeeded2.isNotEmpty()) {
            LOG.error("Irrecoverable fitting error")
            LOG.error("means: " + means.toString())
            LOG.error("ps: " + ps.toString())
            LOG.error("track(s) " + switchNeeded2.joinToString(transform = Int::toString)
                    + " contradict track(s) " + switchNotNeeded2.joinToString(transform = Int::toString))
            throw IllegalStateException("Irrecoverable fitting error")
        }
        if (switchNeeded1.isNotEmpty()) {
            // flip LOW(0) <-> DECREASING(2) and HIGH(3) <-> INCREASING(1)
            for (i in 0 until tracks1) {
                val tmp = model[i + 1]
                model[i + 1] = model[i + tracks1 + 1] as NegBinEmissionScheme
                model[i + tracks1 + 1] = tmp as NegBinEmissionScheme
            }
            probabilityFlip(model, 0, 2, 5)
            probabilityFlip(model, 1, 3, 5)
        }
        if (switchNeeded2.isNotEmpty()) {
            // flip LOW(0) <-> INCREASING(1) and HIGH(3) <-> DECREASING(2)
            for (i in 2 * tracks1 until 2 * tracks1 + tracks2) {
                val tmp = model[i + 1]
                model[i + 1] = model[i + tracks2 + 1] as NegBinEmissionScheme
                model[i + tracks2 + 1] = tmp as NegBinEmissionScheme
            }
            probabilityFlip(model, 0, 1, 5)
            probabilityFlip(model, 2, 3, 5)
        }
        model.updateTransients()
    }


    private fun semanticCheck(fitter: Fitter<MLConstrainedNBHMM>, tracks: Int): Fitter<MLConstrainedNBHMM> {
        return object : Fitter<MLConstrainedNBHMM> by fitter {

            override fun fit(preprocessed: Preprocessed<DataFrame>,
                             title: String,
                             threshold: Double,
                             maxIter: Int): MLConstrainedNBHMM {
                val model = fitter.fit(preprocessed, title, threshold, maxIter)
                stateFlip(model, tracks)
                return model
            }

            override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                             title: String,
                             threshold: Double,
                             maxIter: Int): MLConstrainedNBHMM {
                val model = fitter.fit(preprocessed, title, threshold, maxIter)
                stateFlip(model, tracks)
                return model
            }
        }
    }

    /**
     * Flip states in case when states with HIGH get lower mean than LOW
     */
    @VisibleForTesting
    internal fun stateFlip(model: MLConstrainedNBHMM, tracks: Int) {
        val means = model.means
        val ps = model.successProbabilities
        val switchNeeded = (0 until tracks).filter { means[it] > means[it + tracks] && ps[it] > ps[it + tracks] }
        val switchNotNeeded = (0 until tracks).filter { means[it] < means[it + tracks] || ps[it] < ps[it + tracks] }
        if (switchNeeded.isNotEmpty() && switchNotNeeded.isNotEmpty()) {
            LOG.error("Irrecoverable fitting error")
            LOG.error("means: " + means.toString())
            LOG.error("ps: " + ps.toString())
            LOG.error("track(s) " + switchNeeded.joinToString(transform = Int::toString)
                    + " contradict track(s) " + switchNotNeeded.joinToString(transform = Int::toString))
            throw IllegalStateException("Irrecoverable fitting error")
        }
        if (switchNeeded.isNotEmpty()) {
            // flip LOW(1) <-> HIGH(2)
            for (i in 0 until tracks) {
                val tmp = model[i + 1]
                model[i + 1] = model[i + tracks + 1] as NegBinEmissionScheme
                model[i + tracks + 1] = tmp as NegBinEmissionScheme
            }
            probabilityFlip(model, 1, 2, 3)
        }
        model.updateTransients()
    }

    private fun semanticCheck(fitter: Fitter<MLFreeNBHMM>): Fitter<MLFreeNBHMM> {
        return object : Fitter<MLFreeNBHMM> by fitter {

            override fun fit(preprocessed: Preprocessed<DataFrame>,
                             title: String,
                             threshold: Double,
                             maxIter: Int): MLFreeNBHMM {
                val model = fitter.fit(preprocessed, title, threshold, maxIter)
                stateFlip(model)
                return model
            }

            override fun fit(preprocessed: List<Preprocessed<DataFrame>>,
                             title: String,
                             threshold: Double,
                             maxIter: Int): MLFreeNBHMM {
                val model = fitter.fit(preprocessed, title, threshold, maxIter)
                stateFlip(model)
                return model
            }
        }
    }

    /**
     * Flip states in case when states with HIGH get lower mean than LOW
     */
    @VisibleForTesting
    internal fun stateFlip(model: MLFreeNBHMM) {
        val lowScheme = model[1] as NegBinEmissionScheme
        val highScheme = model[2] as NegBinEmissionScheme
        val meanLow = model.means[0]
        val meanHigh = model.means[1]
        val pLow = lowScheme.successProbability
        val pHigh = highScheme.successProbability
        val meanFlipped = meanLow > meanHigh
        val pFlipped = pLow > pHigh
        if (meanFlipped) {
            LOG.warn("After fitting the model, mean emission in LOW state ($meanLow) is higher than " +
                    "mean emission in HIGH state ($meanHigh).")
        }
        if (pFlipped) {
            LOG.warn("After fitting the model, emission's parameter p in LOW state ($pLow) is higher than " +
                    "emission's parameter p in HIGH state ($pHigh).")
        }
        if (meanFlipped && pFlipped) {
            LOG.warn("This usually indicates that the states were flipped during fitting. We will now flip them back.")
            model[2] = lowScheme
            model[1] = highScheme
            probabilityFlip(model, 1, 2)
        } else if (meanFlipped || pFlipped) {
            LOG.warn("This is generally harmless, but could indicate low quality of data.")
        }
    }

    /**
     * Rearranges the transition and prior probabilities so that the provided states flip
     */
    private fun probabilityFlip(model: MLConstrainedNBHMM, state1: Int, state2: Int, stateNum: Int) {
        for (i in 0 until stateNum) {
            val tmp = model.logTransitionProbabilities[i, state1]
            model.logTransitionProbabilities[i, state1] = model.logTransitionProbabilities[i, state2]
            model.logTransitionProbabilities[i, state2] = tmp
        }
        for (j in 0 until stateNum) {
            val tmp = model.logTransitionProbabilities[state1, j]
            model.logTransitionProbabilities[state1, j] = model.logTransitionProbabilities[state2, j]
            model.logTransitionProbabilities[state2, j] = tmp
        }
        val tmp = model.logPriorProbabilities[state1]
        model.logPriorProbabilities[state1] = model.logPriorProbabilities[state2]
        model.logPriorProbabilities[state2] = tmp
    }

    private fun probabilityFlip(model: MLFreeNBHMM, state1: Int, state2: Int) {
        for (i in 0..2) {
            val tmp = model.logTransitionProbabilities[i, state1]
            model.logTransitionProbabilities[i, state1] = model.logTransitionProbabilities[i, state2]
            model.logTransitionProbabilities[i, state2] = tmp
        }
        for (j in 0..2) {
            val tmp = model.logTransitionProbabilities[state1, j]
            model.logTransitionProbabilities[state1, j] = model.logTransitionProbabilities[state2, j]
            model.logTransitionProbabilities[state2, j] = tmp
        }
        val tmp = model.logPriorProbabilities[state1]
        model.logPriorProbabilities[state1] = model.logPriorProbabilities[state2]
        model.logPriorProbabilities[state2] = tmp
    }


    @JvmOverloads
    fun computeDifferentialPeaks(genomeQuery: GenomeQuery,
                                 queries1: List<InputQuery<Coverage>>,
                                 queries2: List<InputQuery<Coverage>>,
                                 bin: Int = BIN,
                                 fdr: Double = FDR,
                                 gap: Int = GAP): LocationsMergingList {
        val peaks = computeDifferentialPeaksList(genomeQuery, queries1, queries2, bin, fdr, gap)
        return LocationsMergingList.create(genomeQuery, peaks.map { it.location })
    }

    private fun computeDifferentialPeaksList(genomeQuery: GenomeQuery,
                                             queries1: List<InputQuery<Coverage>>,
                                             queries2: List<InputQuery<Coverage>>,
                                             bin: Int = BIN,
                                             fdr: Double = FDR,
                                             gap: Int = GAP): List<Peak> {
        return getDifferentialPeakCallingExperiment(genomeQuery, queries1, queries2, bin)
                .results.getPeaks(genomeQuery, fdr, gap)
    }

    fun computeDirectedDifferencePeaks(genomeQuery: GenomeQuery,
                                       queries1: List<InputQuery<Coverage>>,
                                       queries2: List<InputQuery<Coverage>>,
                                       bin: Int = BIN,
                                       fdr: Double = FDR,
                                       gap: Int = GAP): Pair<List<Peak>, List<Peak>> {
        val spanExperiment = getDifferentialPeakCallingExperiment(genomeQuery, queries1, queries2, bin)
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            spanExperiment.results.getChromosomePeaks(chromosome, fdr, gap, spanExperiment.getData(chromosome))
        }
        val highLow = arrayListOf<Peak>()
        val lowHigh = arrayListOf<Peak>()
        genomeQuery.get().forEach { chromosome ->
            val states = spanExperiment.getStatesDataFrame(chromosome)
            map[chromosome].forEach {
                if (states.getAsFloat(it.startOffset / bin, ZLHID.D.name) >
                        states.getAsFloat(it.startOffset / bin, ZLHID.I.name)) {
                    highLow.add(it)
                } else {
                    lowHigh.add(it)
                }
            }
        }
        return highLow to lowHigh
    }


    fun computeDirectedDifference(genomeQuery: GenomeQuery,
                                  queries1: List<InputQuery<Coverage>>,
                                  queries2: List<InputQuery<Coverage>>,
                                  bin: Int = BIN,
                                  fdr: Double = FDR,
                                  gap: Int = GAP): Pair<LocationsMergingList, LocationsMergingList> {
        val (high, low) = computeDirectedDifferencePeaks(genomeQuery, queries1, queries2, bin, fdr, gap)
        return LocationsMergingList.create(genomeQuery, high.map { it.location }) to
                LocationsMergingList.create(genomeQuery, low.map { it.location })
    }

}