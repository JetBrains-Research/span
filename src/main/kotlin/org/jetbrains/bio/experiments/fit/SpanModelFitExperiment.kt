package org.jetbrains.bio.experiments.fit

import com.google.common.annotations.VisibleForTesting
import com.google.common.math.IntMath
import com.google.gson.GsonBuilder
import kotlinx.support.jdk7.use
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment.Companion.createEffectiveQueries
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.span.CoverageScoresQuery
import org.jetbrains.bio.span.scoresDataFrame
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.statistics.state.ZLHID
import org.jetbrains.bio.util.*
import java.math.RoundingMode
import java.nio.file.Path
import java.util.*
import kotlin.collections.LinkedHashMap

/**
 * Since all the chromosomes are squashed in [SpanModelFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 */
data class SpanFitInformation(val build: String,
                              val data: List<TC>,
                              val labels: List<String>,
                              val fragment: Int?,
                              val binSize: Int,
                              val chromosomesSizes: LinkedHashMap<String, Int>,
                              val version: Int) {

    constructor(genomeQuery: GenomeQuery,
                paths: List<Pair<Path, Path?>>,
                labels: List<String>,
                fragment: Int?,
                binSize: Int) :
            this(genomeQuery.build,
                    paths.map { TC(it.first.toString(), it.second?.toString()) },
                    labels, fragment, binSize,
                    LinkedHashMap<String, Int>().apply {
                        genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
                    },
                    VERSION)

    internal fun checkBuild(build: String) {
        check(this.build == build) {
            "Wrong genome build, expected: ${this.build}, got: $build"
        }
    }

    private fun checkChromosome(chromosome: Chromosome) {
        check(chromosome.name in chromosomesSizes) {
            "Missing chromosome in ${chromosomesSizes.keys.toList()}: ${chromosome.name}"
        }
        check(chromosome.length == chromosomesSizes[chromosome.name]) {
            "Wrong chromosome ${chromosome.name} size, expected ${chromosomesSizes[chromosome.name]}, got: ${chromosome.length}"
        }
    }

    private fun offsetsMap(): IntArray =
            (listOf(0) + chromosomesSizes.keys.sorted().map {
                IntMath.divide(chromosomesSizes[it]!!, binSize, RoundingMode.CEILING)
            }).toIntArray().let {
                Arrays.parallelPrefix(it) { a, b -> a + b }; it
            }


    /**
     * Creates binned offsets for [chromosome] using [binSize]
     */
    fun offsets(chromosome: Chromosome): IntArray {
        checkBuild(chromosome.genome.build)
        checkChromosome(chromosome)
        return chromosome.range.slice(binSize).mapToInt { it.startOffset }.toArray()
    }

    /**
     * Since all chromosomes are squashed into a single data frame for [SpanModelFitExperiment]
     * This method computes indices of data frame, for given [chromosome]
     * See also: [merge] and [split]
     */
    internal fun getChromosomesIndices(chromosome: Chromosome): Pair<Int, Int> {
        checkChromosome(chromosome)
        val offsetsMap = offsetsMap()
        val index = chromosomesSizes.keys.sorted().indexOf(chromosome.name)
        return offsetsMap[index] to offsetsMap[index + 1]
    }

    internal fun save(path: Path) {
        path.parent.createDirectories()
        path.bufferedWriter().use { GSON.toJson(this, it) }
    }

    internal fun merge(statesDataFrame: Map<String, DataFrame>): DataFrame {
        return DataFrame.rowBind(chromosomesSizes.keys.sorted().map { statesDataFrame[it]!! }.toTypedArray())
    }

    internal fun split(dataFrame: DataFrame, genomeQuery: GenomeQuery): Map<String, DataFrame> {
        checkBuild(genomeQuery.build)
        return genomeQuery.get()
                .filter { it.name in chromosomesSizes }
                .map { chromosome ->
                    val (start, end) = getChromosomesIndices(chromosome)
                    chromosome.name to dataFrame.iloc[start until end]
                }.toMap()
    }

    companion object {
        const val VERSION: Int = 1

        /**
         * Using Treatment and Control class instead of [Pair] here for nice GSON serialization
         */
        data class TC(val path: String, val control: String?)


        private val GSON = GsonBuilder()
                .setPrettyPrinting()
                .setFieldNamingStrategy(GSONUtil.NO_MY_UNDESCORE_NAMING_STRATEGY)
                .create()

        fun load(path: Path): SpanFitInformation {
            return path.bufferedReader().use {
                val info = GSON.fromJson(it, SpanFitInformation::class.java)
                checkNotNull(info) {
                    "Failed to load info from $path"
                }
                check(VERSION == info.version) {
                    "Wrong version: expected: $VERSION, got: ${info.version}"
                }
                return@use info
            }
        }
    }
}

data class SpanFitResults(val fitInfo: SpanFitInformation,
                          val model: ClassificationModel,
                          val logNullMemberships: Map<String, DataFrame>)


/**
 * A generic class for [SPAN] (Semi-supervised Peak Analyzer) - tool for analyzing and comparing ChIP-Seq data.
 * Both procedures rely on the Zero Inflated Negative Binomial Restricted Algorithm.
 *
 * It is implemented as [ModelFitExperiment] with different [ClassificationModel] models.
 *
 * Enrichment
 * - States: [ZLH]
 * - Single replicate: [MLFreeNBHMM] zero-inflated HMM with univariate Negative Binomial emissions
 * - Multi replicates: [MLConstrainedNBHMM] zero-inflated HMM with multidimensional Negative Binomial emissions
 *
 * Difference
 * - States: [ZLHID]
 * - Any number of replicates: [MLConstrainedNBHMM]
 */
abstract class SpanModelFitExperiment<out Model : ClassificationModel, State : Any>(
        /** XXX may contain chromosomes without reads, use [genomeQuery] instead. Details: [createEffectiveQueries] */
        externalGenomeQuery: GenomeQuery,
        paths: List<Pair<Path, Path?>>,
        labels: List<String>,
        fragment: Int?,
        val binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        availableStates: Array<State>,
        private val nullHypothesis: NullHypothesis<State>)
    : ModelFitExperiment<Model, State>(
        createEffectiveQueries(externalGenomeQuery, paths, labels, fragment, binSize),
        modelFitter, modelClass, availableStates) {

    val results: SpanFitResults by lazy {
        getOrLoadResults()
    }

    override fun getStatesDataFrame(chromosome: Chromosome): DataFrame = sliceStatesDataFrame(statesDataFrame, chromosome)

    override fun doCalculations() {
        results.logNullMemberships
    }

    val fitInformation = SpanFitInformation(genomeQuery, paths, labels, fragment, binSize)

    private val preprocessedData: List<Preprocessed<DataFrame>> by lazy {
        genomeQuery.get().sortedBy { it.name }.map {
            Preprocessed.of(dataQuery.apply(it))
        }
    }

    // XXX It is important to use get() here, because id is overridden in superclasses
    private val modelPath: Path
        get() = experimentPath / "$id.span"

    private fun calculateModel(): Model {
        MultitaskProgress.addTask(dataQuery.id, (Fitter.MAX_ITERATIONS / 4).toLong())
        val model = modelFitter.fit(preprocessedData, title = dataQuery.id)
        MultitaskProgress.finishTask(dataQuery.id)
        return model
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
                df = df.with("state", model.predict(preprocessed)
                        .map { availableStates[it].toString() }.toTypedArray())
                df
            }.toTypedArray())

    private val statesDataFrame: DataFrame by lazy {
        @Suppress("UNCHECKED_CAST")
        calculateStatesDataFrame(results.model as Model)
    }

    private fun sliceStatesDataFrame(statesDataFrame: DataFrame, chromosome: Chromosome): DataFrame {
        val (start, end) = fitInformation.getChromosomesIndices(chromosome)
        return statesDataFrame.iloc[start until end]
    }

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
                val model = calculateModel()
                model.save(modelPath)

                val informationPath = dir / INFORMATION_JSON
                fitInformation.save(informationPath)

                val statesDataFrame = calculateStatesDataFrame(model)
                val chromosomeToDataFrameMap = genomeQuery.get().associate {
                    val logMemberships = getLogMemberships(sliceStatesDataFrame(statesDataFrame, it))
                    val logNullMemberships = nullHypothesis.apply(logMemberships)
                    // Convert [Double] to [Float] to save space, see #1163
                    it.name to DataFrame().with(NULL, logNullMemberships.toFloatArray())
                }
                val logNullMembershipsDF = fitInformation.merge(chromosomeToDataFrameMap)
                val nullHypothesisPath = dir / NULL_NPZ
                logNullMembershipsDF.save(nullHypothesisPath)

                // Sanity check: model load
                ClassificationModel.load<Model>(modelPath)
                // Sanity check: information load
                SpanFitInformation.load(informationPath)
                val logNullMembershipsMap = fitInformation.split(logNullMembershipsDF, genomeQuery)
                computedResults = SpanFitResults(fitInformation, model, logNullMembershipsMap)
                Tar.compress(p, modelPath.toFile(), informationPath.toFile(), nullHypothesisPath.toFile())
            }
        }
        return if (computedResults != null) {
            LOG.info("Model saved: $modelPath")
            computedResults!!
        } else {
            loadResults(genomeQuery, modelPath).apply {
                check(this.fitInfo.binSize == binSize) {
                    "Wrong bin size: expected $binSize, but got ${this.fitInfo.binSize}"
                }
            }
        }
    }

    companion object {
        private const val INFORMATION_JSON = "information.json"
        private const val MODEL_JSON = "model.json"
        private const val NULL_NPZ = "null.npz"
        const val NULL = "null"


        @VisibleForTesting
        /**
         * Create pair of
         * 1. Effective genomeQuery, i.e. only chromosomes with some reads on them
         * 2. Data query required for [ModelFitExperiment]
         */
        internal fun createEffectiveQueries(genomeQuery: GenomeQuery,
                                            paths: List<Pair<Path, Path?>>,
                                            labels: List<String>,
                                            fragment: Int?,
                                            binSize: Int): Pair<GenomeQuery, Query<Chromosome, DataFrame>> {
            val chromosomes = genomeQuery.get()
            val nonEmptyChromosomes = hashSetOf<Chromosome>()
            paths.forEach { (t, _) ->
                val coverage = ReadsQuery(genomeQuery, t, true, fragment).get()
                nonEmptyChromosomes.addAll(chromosomes.filter { coverage.getBothStrandsCoverage(it.range.on(it)) > 0 })
            }
            chromosomes.filter { it !in nonEmptyChromosomes }.forEach {
                LOG.info("${it.name}: no reads detected, ignoring.")
            }
            val effectiveGenomeQuery = genomeQuery.only(nonEmptyChromosomes.toList().map { it.name }.sorted())
            return effectiveGenomeQuery to object : CachingQuery<Chromosome, DataFrame>() {
                val scores = paths.map { CoverageScoresQuery(effectiveGenomeQuery, it.first, it.second, fragment, binSize) }

                override fun getUncached(input: Chromosome): DataFrame {
                    return scores.scoresDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(scores.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
            }
        }


        fun loadResults(genomeQuery: GenomeQuery, tarPath: Path): SpanFitResults {
            LOG.info("Loading model: $tarPath")
            return withTempDirectory(tarPath.stem) { dir ->
                LOG.debug("Started model file decompress: $tarPath")
                Tar.decompress(tarPath, dir.toFile())

                LOG.debug("Completed model file decompress and started loading: $tarPath")
                val info = SpanFitInformation.load(dir / SpanModelFitExperiment.INFORMATION_JSON)
                // Sanity check
                info.checkBuild(genomeQuery.build)
                val model = ClassificationModel.load<ClassificationModel>(dir / MODEL_JSON)
                val logNullMembershipsDF = DataFrame.load(dir / SpanModelFitExperiment.NULL_NPZ)
                val logNullMembershipsMap = info.split(logNullMembershipsDF, genomeQuery)

                LOG.info("Completed loading model: $tarPath")
                return@withTempDirectory SpanFitResults(info, model, logNullMembershipsMap)
            }
        }
    }
}