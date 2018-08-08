package org.jetbrains.bio.experiments.fit

import com.google.gson.GsonBuilder
import kotlinx.support.jdk7.use
import org.jetbrains.bio.coverage.GenomeScoresQuery
import org.jetbrains.bio.coverage.scoresDataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.NullHypothesis
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.data.DataFrame
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.statistics.state.ZLHID
import org.jetbrains.bio.util.*
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
                              val version: Int = VERSION) {

    constructor(genomeQuery: GenomeQuery,
                paths: List<Pair<Path, Path?>>,
                labels: List<String>,
                fragment: Int?,
                binSize: Int) :
            this(genomeQuery.build,
                    paths.map { TC(it.first.toString(), it.second?.toString()) },
                    labels, fragment, binSize,
                    chromosomesSizes = LinkedHashMap<String, Int>().apply {
                        genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
                    })

    fun checkBinSize(binSize: Int) {
        check(this.binSize == binSize) {
            "Wrong bin size, expected: ${this.binSize}, got: $binSize"
        }
    }

    private fun checkBuild(build: String) {
        check(this.build == build) {
            "Wrong genome build, expected: ${this.build}, got: $build"
        }
    }

    fun checkGenomeQuery(genomeQuery: GenomeQuery) {
        checkBuild(genomeQuery.build)
        val names = genomeQuery.get().map { it.name }.sorted()
        check(chromosomesSizes.keys.toList() == names) {
            "Wrong chromosomes, expected: ${chromosomesSizes.keys.toList()}, got: $names"
        }
        genomeQuery.get().forEach { checkChromosome(it) }
    }

    private fun checkChromosome(chromosome: Chromosome) {
        check(chromosome.name in chromosomesSizes) {
            "Missing chromosome in ${chromosomesSizes.keys.toList()}: ${chromosome.name}"
        }
        check(chromosome.length == chromosomesSizes[chromosome.name]) {
            "Wrong chromosome ${chromosome.name} size, expected ${chromosomesSizes[chromosome.name]}, got: ${chromosome.length}"
        }
    }

    fun sortedChromosomes(genomeQuery: GenomeQuery): List<Chromosome> {
        checkGenomeQuery(genomeQuery)
        return genomeQuery.get().sortedBy { it.name }
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
     */
    internal fun getChromosomesIndices(genomeQuery: GenomeQuery, chromosome: Chromosome): Pair<Int, Int> {
        checkGenomeQuery(genomeQuery)
        val offsets = (listOf(0) + sortedChromosomes(genomeQuery).map { offsets(it).size })
                .toIntArray().let {
                    Arrays.parallelPrefix(it) { a, b -> a + b }; it
                }
        checkGenomeQuery(genomeQuery)
        val sortedChromosomes = sortedChromosomes(genomeQuery)
        val index = (0 until sortedChromosomes.size).find {
            sortedChromosomes[it] == chromosome
        }
        check(index != null) {
            "Failed to find chromosome ${chromosome.name} in $sortedChromosomes"
        }
        return offsets[index!!] to offsets[index + 1]
    }

    fun save(path: Path) {
        path.parent.createDirectories()
        path.bufferedWriter().use { GSON.toJson(this, it) }
    }

    fun dataFrameToMap(dataFrame: DataFrame, genomeQuery: GenomeQuery): Map<String, DataFrame> {
        checkGenomeQuery(genomeQuery)
        return genomeQuery.get().map { chromosome ->
            val (start, end) = getChromosomesIndices(genomeQuery, chromosome)
            chromosome.name to dataFrame.iloc[start until end]
        }.toMap()
    }

    companion object {
        const val VERSION: Int = 1

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
                    "Wrong $path version: expected: $VERSION, got: ${info.version}"
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
        genomeQuery: GenomeQuery,
        paths: List<Pair<Path, Path?>>,
        labels: List<String>,
        fragment: Int?,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        availableStates: Array<State>,
        private val nullHypothesis: NullHypothesis<State>)
    : ModelFitExperiment<Model, State>(
        genomeQuery,
        createDataQuery(genomeQuery, paths, labels, fragment, binSize),
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
        fitInformation.sortedChromosomes(genomeQuery).map {
            Preprocessed.of(dataQuery.apply(it))
        }
    }

    // XXX It is important to use get() here, because id is overridden in superclasses
    private val tarPath: Path
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
        val (start, end) = fitInformation.getChromosomesIndices(genomeQuery, chromosome)
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
        tarPath.checkOrRecalculate("Model fit") { (p) ->
            withTempDirectory(tarPath.stem) { dir ->
                val modelPath = dir / MODEL_JSON
                val model = calculateModel()
                model.save(modelPath)

                val informationPath = dir / INFORMATION_JSON
                fitInformation.save(informationPath)

                val statesDataFrame = calculateStatesDataFrame(model)

                val nullHypothesisPath = dir / NULL_NPZ
                val logNullMembershipsDF = DataFrame.rowBind(
                        fitInformation.sortedChromosomes(genomeQuery).map { chromosome ->
                            val logMemberships = getLogMemberships(sliceStatesDataFrame(statesDataFrame, chromosome))
                            val logNullMemberships = nullHypothesis.apply(logMemberships)
                            // Convert [Double] to [Float] to save space, see #1163
                            DataFrame().with(NULL, logNullMemberships.toFloatArray())
                        }.toTypedArray()
                )
                logNullMembershipsDF.save(nullHypothesisPath)

                // Sanity check: model load
                ClassificationModel.load<Model>(modelPath)
                // Sanity check: information load
                SpanFitInformation.load(informationPath)
                val logNullMembershipsMap = fitInformation.dataFrameToMap(logNullMembershipsDF, genomeQuery)
                computedResults = SpanFitResults(fitInformation, model, logNullMembershipsMap)
                Tar.compress(p, modelPath.toFile(), informationPath.toFile(), nullHypothesisPath.toFile())
            }
        }
        return if (computedResults != null) {
            LOG.info("Model saved: $tarPath")
            computedResults!!
        } else {
            loadResults(genomeQuery, tarPath)
        }
    }

    companion object {
        private const val INFORMATION_JSON = "information.json"
        private const val MODEL_JSON = "model.json"
        private const val NULL_NPZ = "null.npz"
        const val NULL = "null"


        private fun createDataQuery(genomeQuery: GenomeQuery,
                                    paths: List<Pair<Path, Path?>>,
                                    labels: List<String>,
                                    fragment: Int?,
                                    binSize: Int): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                val scores = paths.map { GenomeScoresQuery(genomeQuery, it.first, it.second, fragment, binSize) }

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
                info.checkGenomeQuery(genomeQuery)
                val model = ClassificationModel.load<ClassificationModel>(dir / MODEL_JSON)
                val logNullMembershipsDF = DataFrame.load(dir / SpanModelFitExperiment.NULL_NPZ)
                val logNullMembershipsMap = info.dataFrameToMap(logNullMembershipsDF, genomeQuery)

                LOG.info("Completed loading model: $tarPath")
                return@withTempDirectory SpanFitResults(info, model, logNullMembershipsMap)
            }
        }
    }
}