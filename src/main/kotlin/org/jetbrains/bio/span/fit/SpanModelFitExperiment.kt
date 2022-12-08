package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiment.Experiment
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.ReadsQuery
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
 * @param fixedModelPath    If not null, the experiment will use this path for saving/loading the results.
 *                          Otherwise, [defaultModelPath] will be used (it usually depends on [fitInformation] id).
 *
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
    private val maxIterations: Int = Fitter.MAX_ITERATIONS,
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
        return modelFitter.fit(
            preprocessedData,
            title = dataQuery.id,
            threshold = threshold,
            maxIterations = maxIterations
        )
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
                LOG.info("Computing data model...")
                val model = calculateModel()
                LOG.info("Done computing data model")
                LOG.info("Saving model information")
                model.save(modelPath)
                LOG.debug("Model saved to $modelPath")

                val informationPath = dir / INFORMATION_JSON
                fitInformation.save(informationPath)
                LOG.debug("Fit information saved to $informationPath")

                LOG.info("Analyzing model bins enrichment")
                val statesDataFrame = calculateStatesDataFrame(model)
                val chromosomes = genomeQuery.get()
                val chromosomeToDataFrameMap = chromosomes.associate {
                    val logMemberships = getLogMemberships(sliceStatesDataFrame(statesDataFrame, it))
                    val logNullMemberships = nullHypothesis.apply(logMemberships)
                    // Convert [Double] to [Float] to save space, see #1163
                    it.name to DataFrame().with(NULL, logNullMemberships.toFloatArray())
                }
                val logNullMembershipsDF = fitInformation.merge(chromosomeToDataFrameMap)
                LOG.debug("Done null hypothesis log memberships")

                val logNullMembershipsPath = dir / NULL_NPZ
                logNullMembershipsDF.save(logNullMembershipsPath)
                LOG.debug("LogNullMemberships saved to $logNullMembershipsPath")
                var statesDataFrameMap: Map<String, DataFrame>? = null
                if (saveExtendedInfo) {
                    LOG.debug("Saving full states dataframe")
                    val statesDataFramePath = dir / "states.npz"
                    statesDataFrame.save(statesDataFramePath)
                    LOG.debug("States saved to $statesDataFramePath")
                    Tar.compress(
                        p,
                        *listOf(
                            modelPath,
                            informationPath,
                            logNullMembershipsPath,
                            statesDataFramePath
                        ).map(Path::toFile)
                            .toTypedArray()
                    )
                    statesDataFrameMap = fitInformation.split(statesDataFrame, genomeQuery)
                } else {
                    Tar.compress(
                        p,
                        *listOf(modelPath, informationPath, logNullMembershipsPath).map(Path::toFile)
                            .toTypedArray()
                    )
                }

                val logNullMembershipsMap = fitInformation.split(logNullMembershipsDF, genomeQuery)
                computedResults = SpanFitResults(fitInformation, model, logNullMembershipsMap, statesDataFrameMap)
                LOG.info("Done saving model.")

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
                    ReadsQuery(genomeQuery, c, unique, fragment, showLibraryInfo = false).get()
                }
                nonEmptyChromosomes.addAll(chromosomes.filter { coverage.getBothStrandsCoverage(it.range.on(it)) > 0 })
            }

            if (nonEmptyChromosomes.isEmpty()) {
                val errMessage = "Model can't be trained on empty coverage, exiting."
                LOG.error(errMessage)
                throw IllegalStateException(errMessage)
            }

            val emptyChromosomes = chromosomes.filter { it !in nonEmptyChromosomes }
            if (emptyChromosomes.isNotEmpty()) {
                LOG.info("Chromosomes with no reads detected are ignored. Use --debug for details.")
                LOG.debug("Ignored chromosomes: ${emptyChromosomes.joinToString(",") { it.name }}")
            }

            return GenomeQuery(genomeQuery.genome, *nonEmptyChromosomes.map { it.name }.toTypedArray())
        }


        fun loadResults(
            genomeQuery: GenomeQuery? = null,
            tarPath: Path
        ): SpanFitResults {
            LOG.info("Loading model: $tarPath")
            return withTempDirectory(tarPath.stem) { dir ->
                LOG.debug("Started model file decompress: ${tarPath.stem}")
                Tar.decompress(tarPath, dir.toFile())

                LOG.debug("Completed model file decompress and started loading: ${tarPath.stem}")
                val info = SpanFitInformation.load<SpanFitInformation>(dir / INFORMATION_JSON)
                // Check genome build
                genomeQuery?.let { info.checkGenome(it.genome) }
                // Load model and PEPs
                val model = ClassificationModel.load<ClassificationModel>(dir / MODEL_JSON)
                val logNullMembershipsDF = DataFrame.load(dir / NULL_NPZ)
                val logNullMembershipsMap = info.split(logNullMembershipsDF, genomeQuery)
                var statesDfMap: Map<String, DataFrame>? = null
                val statesPath = dir / "states.npz"
                LOG.info("Loading states data frame ${tarPath.stem}")
                if (statesPath.exists) {
                    statesDfMap = info.split(DataFrame.load(statesPath), genomeQuery)
                }
                LOG.info("Completed loading model: ${tarPath.stem}")
                return@withTempDirectory SpanFitResults(info, model, logNullMembershipsMap, statesDfMap)
            }
        }
    }
}

