package org.jetbrains.bio.experiments.fit

import com.google.gson.GsonBuilder
import kotlinx.support.jdk7.use
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.*
import org.jetbrains.bio.statistics.*
import org.jetbrains.bio.statistics.data.DataFrame
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.util.*
import java.nio.file.Path
import java.util.*
import kotlin.collections.LinkedHashMap


/**
 * Since all the chromosomes are squashed in [CoverageFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 */
data class CoverageFitInformation(val description: String,
                                  val binSize: Int,
                                  val build: String,
                                  val chromosomesSizes: LinkedHashMap<String, Int>) {
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
     * Since all chromosomes are squashed into a single data frame for [CoverageFitExperiment]
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

        fun of(description: String, binSize: Int, genomeQuery: GenomeQuery): CoverageFitInformation {
            return CoverageFitInformation(
                    description = description,
                    binSize = binSize,
                    build = genomeQuery.build,
                    chromosomesSizes = LinkedHashMap<String, Int>().apply {
                        genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
                    })
        }

        private val GSON = GsonBuilder()
                .setPrettyPrinting()
                .setFieldNamingStrategy(GSONUtil.NO_MY_UNDESCORE_NAMING_STRATEGY)
                .create()

        fun load(path: Path): CoverageFitInformation {
            return path.bufferedReader().use {
                val info = GSON.fromJson(it, CoverageFitInformation::class.java)
                checkNotNull(info) {
                    "failed to load info from $path"
                }
                return@use info
            }
        }
    }
}

data class CoverageFitResults(val fitInfo: CoverageFitInformation,
                              val model: ClassificationModel,
                              val logNullMemberships: Map<String, DataFrame>)

/**
 * Creates simple [DataFrame] for binned coverage given [chromosome] and [binSize].
 * See call site for exact [labels].
 */
internal fun List<InputQuery<Coverage>>.coverageDataFrame(chromosome: Chromosome,
                                                          binSize: Int,
                                                          labels: Array<String>): DataFrame {
    var res = DataFrame()
    forEachIndexed { d, inputQuery ->
        val coverage = inputQuery.get()
        val binnedCoverage = coverage.getBinnedChromosomeCoverage(chromosome, binSize).toIntArray()
        res = res.with(labels[d], binnedCoverage)
    }
    return res
}

/**
 * A generic experiment for evaluating binned [Coverage] classification models.
 * All the chromosomes are processed by single model.
 */
abstract class CoverageFitExperiment<out Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        private val coverageQuery: Query<Chromosome, DataFrame>,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        availableStates: Array<State>,
        private val nullHypothesis: NullHypothesis<State>)
    : ModelFitExperiment<Model, State>(genomeQuery, coverageQuery, modelFitter, modelClass, availableStates) {

    val results: CoverageFitResults by lazy {
        getOrLoadResults()
    }

    override fun getStatesDataFrame(chromosome: Chromosome): DataFrame = sliceStatesDataFrame(statesDataFrame, chromosome)

    override fun doCalculations() {
        results.logNullMemberships
    }

    private val fitInformation = CoverageFitInformation.of(coverageQuery.description, binSize, genomeQuery)

    private val preprocessedData: List<Preprocessed<DataFrame>> by lazy {
        fitInformation.sortedChromosomes(genomeQuery).map {
            Preprocessed.of(coverageQuery.apply(it))
        }
    }

    @Suppress("LeakingThis")
    private val tarPath: Path = experimentPath / "$id.span"

    private fun calculateModel(): Model {
        MultitaskProgress.addTask(coverageQuery.id, (Fitter.MAX_ITERATIONS / 4).toLong())
        val model = modelFitter.fit(preprocessedData, title = coverageQuery.id)
        MultitaskProgress.finishTask(coverageQuery.id)
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
     * Compute and save [CoverageFitResults], i.e. fit information, trained model and null hypothesis probabilities.
     * If already processed, load them [loadResults].
     *
     * IMPORTANT!
     * We take care not to access any of the lazy properties here, since they depend on this method for initialization.
     */
    private fun getOrLoadResults(): CoverageFitResults {
        var computedResults: CoverageFitResults? = null
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
                        fitInformation.sortedChromosomes(this.genomeQuery).map { chromosome ->
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
                CoverageFitInformation.load(informationPath)
                val logNullMembershipsMap = fitInformation.dataFrameToMap(logNullMembershipsDF, genomeQuery)
                computedResults = CoverageFitResults(fitInformation, model, logNullMembershipsMap)
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

        fun loadResults(genomeQuery: GenomeQuery, tarPath: Path): CoverageFitResults {
            LOG.info("Loading model: $tarPath")
            return withTempDirectory(tarPath.stem) { dir ->
                LOG.debug("Started model file decompress: $tarPath")
                Tar.decompress(tarPath, dir.toFile())

                LOG.debug("Completed model file decompress and started loading: $tarPath")
                val info = CoverageFitInformation.load(dir / CoverageFitExperiment.INFORMATION_JSON)
                // Sanity check
                info.checkGenomeQuery(genomeQuery)
                val model = ClassificationModel.load<ClassificationModel>(dir / MODEL_JSON)
                val logNullMembershipsDF = DataFrame.load(dir / CoverageFitExperiment.NULL_NPZ)
                val logNullMembershipsMap = info.dataFrameToMap(logNullMembershipsDF, genomeQuery)

                LOG.info("Completed loading model: $tarPath")
                return@withTempDirectory CoverageFitResults(info, model, logNullMembershipsMap)
            }
        }
    }
}

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanPeakCallingExperiment<Model : ClassificationModel, State : Any> : CoverageFitExperiment<Model, State> {

    constructor(genomeQuery: GenomeQuery,
                coverageQuery: InputQuery<Coverage>,
                modelFitter: Fitter<Model>,
                modelClass: Class<Model>,
                binSize: Int,
                states: Array<State>,
                nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery, createDataQuery(binSize, coverageQuery),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    constructor(genomeQuery: GenomeQuery,
                coverageQueries: List<InputQuery<Coverage>>,
                modelFitter: Fitter<Model>,
                modelClass: Class<Model>,
                binSize: Int,
                states: Array<State>,
                nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery,
                    createDataQuery(binSize, coverageQueries),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    override val id: String
        get() = dataQuery.id

    companion object {
        const val X_PREFIX = "x"
        const val D_PREFIX = "d"

        fun guessReadsQueries(genomeQuery: GenomeQuery, description: String): List<InputQuery<Coverage>> {
            return description.substringBefore(" analysis bin size:").split(";").mapNotNull {
                ReadsQuery.guessByDescription(genomeQuery, it)
            }
        }

        private fun createDataQuery(binSize: Int,
                                    coverageQuery: InputQuery<Coverage>): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return listOf(coverageQuery).coverageDataFrame(input, binSize, arrayOf(X_PREFIX))
                }

                override val id: String get() = reduceIds(listOf(coverageQuery.id, "$binSize"))

                override val description: String get() = "${coverageQuery.description}; analysis bin size: $binSize"
            }
        }

        private fun createDataQuery(binSize: Int,
                                    coverageQueries: List<InputQuery<Coverage>>): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return coverageQueries.coverageDataFrame(input, binSize,
                            MultiLabels.generate(D_PREFIX, coverageQueries.size))
                }

                override val id: String
                    get() = reduceIds(coverageQueries.map { it.id } + listOf("$binSize"))

                override val description: String
                    get() = "${coverageQueries.joinToString(";") { it.description }} analysis bin size: $binSize"
            }
        }
    }
}

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment<Model : ClassificationModel, State : Any> : CoverageFitExperiment<Model, State> {

    constructor(genomeQuery: GenomeQuery,
                coverageQueries1: List<InputQuery<Coverage>>,
                coverageQueries2: List<InputQuery<Coverage>>,
                binSize: Int, modelFitter: Fitter<Model>, modelClass: Class<Model>,
                states: Array<State>, nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery, createDataQuery(binSize, coverageQueries1, coverageQueries2),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    constructor(genomeQuery: GenomeQuery,
                coverageQuery1: InputQuery<Coverage>,
                coverageQuery2: InputQuery<Coverage>,
                binSize: Int,
                modelFitter: Fitter<Model>,
                modelClass: Class<Model>,
                states: Array<State>, nullHypothesis: NullHypothesis<State>) :
            super(genomeQuery, createDataQuery(binSize, coverageQuery1, coverageQuery2),
                    binSize, modelFitter, modelClass, states, nullHypothesis)

    override val id: String
        get() = "${dataQuery.id}_diff"

    companion object {
        const val TRACK1_PREFIX = "track1_"
        const val TRACK2_PREFIX = "track2_"

        internal fun createDataQuery(binSize: Int,
                                     coverageQueries1: List<InputQuery<Coverage>>,
                                     coverageQueries2: List<InputQuery<Coverage>>): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    val labels1 = MultiLabels.generate(TRACK1_PREFIX, coverageQueries1.size)
                    val res1 = coverageQueries1.coverageDataFrame(input, binSize, labels1)
                    val labels2 = MultiLabels.generate(TRACK2_PREFIX, coverageQueries2.size)
                    val res2 = coverageQueries2.coverageDataFrame(input, binSize, labels2)
                    return DataFrame.columnBind(res1, res2)
                }

                override val id: String
                    get() = reduceIds((coverageQueries1 + coverageQueries2).map { it.id } + "$binSize")

                override val description: String
                    get() = "${coverageQueries1.joinToString(", ") { it.description }} " +
                            "compared to ${coverageQueries1.joinToString(", ") { it.description }} bin size: $binSize"
            }
        }

        internal fun createDataQuery(binSize: Int,
                                     coverageQuery1: InputQuery<Coverage>,
                                     coverageQuery2: InputQuery<Coverage>): Query<Chromosome, DataFrame> {
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    val labels = arrayOf(coverageQuery1.id, coverageQuery2.id)
                    return listOf(coverageQuery1, coverageQuery2).coverageDataFrame(input, binSize, labels)
                }

                override val id: String
                    get() = reduceIds(listOf(coverageQuery1, coverageQuery2).map { it.id } + listOf("$binSize"))

                override val description: String
                    get() =
                        ("${coverageQuery1.description} compared to ${coverageQuery2.description} bin size: $binSize")
            }
        }
    }
}