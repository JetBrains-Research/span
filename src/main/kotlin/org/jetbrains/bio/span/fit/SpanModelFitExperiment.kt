package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiment.Experiment
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.span.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.span.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.span.statistics.mixture.PoissonRegressionMixture
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.f64Array
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.jetbrains.bio.statistics.model.Fitter
import org.jetbrains.bio.statistics.toFloatArray
import org.jetbrains.bio.util.*
import org.jetbrains.bio.viktor.F64Array
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path


/**
 * A generic class for Span (Semi-supervised Peak ANalyzer), a tool for analyzing and comparing ChIP-Seq data.
 *
 * Span can utilize various models ([Model]) and inputs ([FitInfo]),
 * as reflected by the class's generic and abstract nature.
 *
 * The end result of the experiment is the [results] property. It's lazy (won't do any calculation until
 * actually accessed) and cached (if possible, will be loaded from the previously created file, if not, will
 * be saved to a file after the computation). The results are saved in a TAR file.
 *
 * Current implementations:
 * - [SpanPeakCallingExperiment] -- enrichment analysis (peak calling).
 *   - States: [ZLH]
 *   - Fit information: [Span1AnalyzeFitInformation]
 *   - Single replicate: [MLFreeNBHMM] zero-inflated HMM with univariate negative binomial emissions
 *   - Multi replicates: [MLConstrainedNBHMM] zero-inflated HMM with multidimensional negative binomial emissions *
 * - [SpanDifferentialPeakCallingExperiment] -- enrichment comparison (differential peak calling).
 *   - States: [ZLHID]
 *   - Fit information: [Span1CompareFitInformation]
 *   - Any number of replicates: [MLConstrainedNBHMM] zero-inflated HMM with multidimensional
 *   negative binomial emissions *
 * - [Span2PeakCallingExperiment] -- enrichment analysis (peak calling).
 *   - States: [ZLH]
 *   - Fit information: [Span2FitInformation]
 *   - Single replicate: [PoissonRegressionMixture] a mixture of Poisson GLMs
 *
 * @param fixedModelPath If not null, the experiment will use this path for saving/loading the results. Otherwise,
 * [defaultModelPath] will be used (it usually depends of [fitInformation] id).
 */
abstract class SpanModelFitExperiment<
        out Model : ClassificationModel, out FitInfo : SpanFitInformation, State : Any
        > protected constructor(
    val fitInformation: FitInfo,
    private val modelFitter: Fitter<Model>,
    private val modelClass: Class<out Model>,
    private val availableStates: Array<State>,
    private val nullHypothesis: NullHypothesis<State>,
    private val fixedModelPath: Path? = null,
    private val threshold: Double = Fitter.THRESHOLD,
    private val maxIter: Int = Fitter.MAX_ITERATIONS,
    private val saveExtendedInfo: Boolean = false
) : Experiment("fit") {

    val genomeQuery = fitInformation.genomeQuery()
    private val dataQuery = fitInformation.dataQuery

    /**
     * Preprocessed data by chromosomes, chromosomes are sorted by name.
     */
    private val preprocessedData: List<Preprocessed<DataFrame>> by lazy {
        genomeQuery.get().sortedBy { it.name }.map { Preprocessed.of(dataQuery.apply(it)) }
    }

    val results: SpanFitResults by lazy {
        getOrLoadResults()
    }

    fun getStatesDataFrame(chromosome: Chromosome): DataFrame = sliceStatesDataFrame(statesDataFrame, chromosome)

    override fun doCalculations() {
        results.logNullMemberships
    }

    /**
     * A unique path used to store/load results if [fixedModelPath] is null.
     *
     * This property is normally implemented through [fitInformation] id.
     */
    abstract val defaultModelPath: Path

    /**
     * We use "get" because we need the actual value of [defaultModelPath] implemented in the descendant class.
     */
    private val modelPath get() = fixedModelPath ?: defaultModelPath

    private fun calculateModel(): Model {
        return modelFitter.fit(preprocessedData, title = dataQuery.id, threshold = threshold, maxIter = maxIter)
    }

    private fun calculateStatesDataFrame(model: Model): DataFrame = DataFrame.rowBind(
        preprocessedData.map { preprocessed ->
            val logMemberships = model.evaluate(preprocessed)
            var df = DataFrame()
            availableStates.forEachIndexed { j, state ->
                val f64Array = logMemberships.V[j]
                // Convert [Double] to [Float] to save space, see #1163
                df = df.with(state.toString(), f64Array.toFloatArray())
            }
            df = df.with(
                "state", model.predict(preprocessed)
                    .map { availableStates[it].toString() }.toTypedArray()
            )
            df
        }.toTypedArray()
    )

    private val statesDataFrame: DataFrame by lazy {
        @Suppress("UNCHECKED_CAST")
        calculateStatesDataFrame(results.model as Model)
    }

    /**
     * The [statesDataFrame] contains aggregated information for the whole genome,
     * this method returns data chunk for given [chromosome].
     * @return dataframe slice for given chromosome.
     */
    private fun sliceStatesDataFrame(statesDataFrame: DataFrame, chromosome: Chromosome): DataFrame {
        val (start, end) = fitInformation.getChromosomesIndices(chromosome)
        return statesDataFrame.iloc[start until end]
    }

    /**
     * Return map state -> f64 array of log probabilities for each position to be in given hidden state.
     * Membership = Probability here.
     */
    private fun getLogMemberships(chromosomeStatesDF: DataFrame): Map<State, F64Array> =
        availableStates.associateBy({ it }) { chromosomeStatesDF.f64Array(it.toString()) }

    /**
     * Compute and save [SpanFitResults], i.e. fit information, trained model and null hypothesis probabilities.
     * If already processed, load them [loadResults].
     *
     * IMPORTANT!
     * We take care not to access any of the lazy properties here, since they depend on this method for initialization.
     */
    private fun getOrLoadResults(): SpanFitResults {
        var computedResults: SpanFitResults? = null
        modelPath.checkOrRecalculate("Model fit") { (p) ->
            withTempDirectory(modelPath.stem) { dir ->
                val modelPath = dir / MODEL_JSON
                LOG.info("Computing data model")
                val model = calculateModel()
                LOG.info("Done computing data model")
                LOG.info("Saving model information...")
                model.save(modelPath)
                LOG.debug("Model saved to $modelPath")
                ClassificationModel.load<Model>(modelPath)
                LOG.debug("Sanity check: model loaded")

                val informationPath = dir / INFORMATION_JSON
                fitInformation.save(informationPath)
                LOG.debug("Fit information saved to $informationPath")
                SpanFitInformation.load<SpanFitInformation>(informationPath)
                LOG.debug("Sanity check: information loaded")

                LOG.debug("Computing states dataframe")
                val statesDataFrame = calculateStatesDataFrame(model)
                LOG.debug("Computing null hypothesis log memberships")
                val chromosomes = genomeQuery.get()
                val chromosomeToDataFrameMap = chromosomes.associate {
                    val logMemberships = getLogMemberships(sliceStatesDataFrame(statesDataFrame, it))
                    val logNullMemberships = nullHypothesis.apply(logMemberships)
                    // Convert [Double] to [Float] to save space, see #1163
                    it.name to DataFrame().with(NULL, logNullMemberships.toFloatArray())
                }
                val logNullMembershipsDF = fitInformation.merge(chromosomeToDataFrameMap)
                LOG.debug("Done null hypothesis log memberships")

                val nullHypothesisPath = dir / NULL_NPZ
                logNullMembershipsDF.save(nullHypothesisPath)
                LOG.debug("LogNullMemberships saved to $nullHypothesisPath")
                val logNullMembershipsMap = fitInformation.split(logNullMembershipsDF, genomeQuery)

                if (saveExtendedInfo) {
                    val extendedPaths = arrayListOf<Path>()
                    LOG.info("Extended information in enabled, preparing additional information.")
                    val coverageDfMap = hashMapOf<String, DataFrame>()
                    preprocessedData.forEachIndexed { i, preprocessed ->
                        val chromosome = chromosomes[i]
                        val preprocessedPath = dir / "coverage_${chromosome.name}.npz"
                        val covDf = preprocessed.get()
                        coverageDfMap[chromosome.name] = covDf
                        covDf.save(preprocessedPath)
                        LOG.debug("Preprocessed coverage ${chromosome.name} saved to $preprocessedPath")
                        extendedPaths.add(preprocessedPath)
                    }
                    LOG.debug("Saving full states dataframe")
                    val statesDataFramePath = dir / "states.npz"
                    statesDataFrame.save(statesDataFramePath)
                    LOG.debug("States saved to $statesDataFramePath")
                    extendedPaths.add(statesDataFramePath)
                    Tar.compress(
                        p, *(listOf(modelPath, informationPath, nullHypothesisPath) + extendedPaths).map(Path::toFile)
                            .toTypedArray()
                    )
                    computedResults = SpanFitResultsExt(
                        fitInformation, model, logNullMembershipsMap,
                        coverageDfMap, fitInformation.split(statesDataFrame, genomeQuery)
                    )
                    LOG.info("Done saving model with additional information.")
                } else {
                    Tar.compress(
                        p, *(listOf(modelPath, informationPath, nullHypothesisPath)).map(Path::toFile)
                            .toTypedArray()
                    )
                    computedResults = SpanFitResults(fitInformation, model, logNullMembershipsMap)
                    LOG.info("Done saving model.")
                }
            }
        }
        return if (computedResults != null) {
            LOG.info("Model saved: $modelPath")
            computedResults!!
        } else {
            val loadedResults = loadResults(genomeQuery, modelPath)
            val loadedFitInfo = loadedResults.fitInfo
            check(loadedFitInfo == fitInformation) {
                "Wrong model information: expected $fitInformation, but got $loadedFitInfo"
            }
            loadedResults
        }
    }

    companion object {
        private const val INFORMATION_JSON = "information.json"
        private const val MODEL_JSON = "model.json"
        private const val NULL_NPZ = "null.npz"
        const val NULL = "null"

        val LOG: Logger = LoggerFactory.getLogger(SpanModelFitExperiment::class.java)

        /**
         * Retain only the chromosomes for which at least one treatment file has at least one read on them.
         */
        fun filterGenomeQueryWithData(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            fragment: Fragment,
            unique: Boolean = true
        ): GenomeQuery {
            val chromosomes = genomeQuery.get()
            val nonEmptyChromosomes = hashSetOf<Chromosome>()
            paths.forEach { (t, c) ->
                val coverage = ReadsQuery(genomeQuery, t, unique, fragment, showLibraryInfo = false).get()
                if (c != null) {
                    // we have to be sure that the control coverage cache is calculated for the full genome query,
                    // otherwise we can get some very hard-to-catch bugs later
                    ReadsQuery(genomeQuery, c, unique, fragment).get()
                }
                nonEmptyChromosomes.addAll(chromosomes.filter { coverage.getBothStrandsCoverage(it.range.on(it)) > 0 })
            }
            chromosomes.filter { it !in nonEmptyChromosomes }.forEach {
                LOG.info("${it.name}: no reads detected, ignoring.")
            }
            if (nonEmptyChromosomes.isEmpty()) {
                val errMessage = "Model can't be trained on empty coverage, exiting."
                LOG.error(errMessage)
                throw IllegalStateException(errMessage)
            }
            return GenomeQuery(genomeQuery.genome, *nonEmptyChromosomes.map { it.name }.toTypedArray())
        }


        fun loadResults(genomeQuery: GenomeQuery? = null, tarPath: Path): SpanFitResults {
            LOG.info("Loading model: $tarPath")
            return withTempDirectory(tarPath.stem) { dir ->
                LOG.debug("Started model file decompress: $tarPath")
                Tar.decompress(tarPath, dir.toFile())

                LOG.debug("Completed model file decompress and started loading: $tarPath")
                val info = SpanFitInformation.load<SpanFitInformation>(dir / INFORMATION_JSON)
                // Check genome build
                genomeQuery?.let { info.checkGenome(it.genome) }
                // Load model and PEPs
                val model = ClassificationModel.load<ClassificationModel>(dir / MODEL_JSON)
                val logNullMembershipsDF = DataFrame.load(dir / NULL_NPZ)
                val logNullMembershipsMap = info.split(logNullMembershipsDF, genomeQuery)

                val statesPath = dir / "states.npz"
                if (statesPath.exists) {
                    LOG.info("Loading extended model file")
                    LOG.info("Loading states data frame")
                    val statesDfMap = info.split(DataFrame.load(statesPath), genomeQuery)
                    val coveragesDfMap = hashMapOf<String, DataFrame>()
                    val chromosomeNames = if (genomeQuery != null)
                        genomeQuery.get().map { it.name }
                    else
                        info.chromosomesSizes.keys
                    chromosomeNames.forEach { chromosome ->
                        val chrCoveragePath = dir / "coverage_$chromosome.npz"
                        if (chrCoveragePath.exists) {
                            LOG.info("Loading coverage data for $chromosome")
                            coveragesDfMap[chromosome] = DataFrame.load(chrCoveragePath)
                        } else {
                            LOG.info("No coverage information available for $chromosome")
                        }
                    }
                    LOG.info("Completed loading extended model: $tarPath")
                    return SpanFitResultsExt(info, model, logNullMembershipsMap, coveragesDfMap, statesDfMap)
                }

                LOG.info("Completed loading model: $tarPath")
                return@withTempDirectory SpanFitResults(info, model, logNullMembershipsMap)
            }
        }
    }
}

