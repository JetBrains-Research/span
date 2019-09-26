package org.jetbrains.bio.experiments.fit

import com.google.common.math.IntMath
import com.google.gson.*
import com.google.gson.reflect.TypeToken
import kotlinx.support.jdk7.use
import org.apache.log4j.Logger
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiment.Experiment
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hmm.MLFreeNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.mixture.PoissonRegressionMixture
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.statistics.state.ZLHID
import org.jetbrains.bio.util.*
import org.jetbrains.bio.viktor.F64Array
import java.lang.reflect.Type
import java.math.RoundingMode
import java.nio.file.Path
import java.util.*

/**
 * The most common interface for all fit information classes.
 *
 * [SpanFitInformation] instance is designed to contain all information necessary to uniquely identify the input
 * of a Span-like model fitting experiment. For example, [Span1AnalyzeFitInformation] completely describes
 * the input of the classical `span analyze` command.
 *
 * [SpanFitInformation] object is a part of [SpanFitResults], and its type is type parameter
 * of [SpanModelFitExperiment].
 *
 * All Span-like experiments produce a single squashed float array of log null probabilities ("null.npz").
 * This interface contains methods to squash ([merge]) and unsquash ([split]) the chromosome-wise dataframes.
 * It can also generate bin start [offsets] for a single chromosome.
 *
 * @property build Genome build (assembly).
 * @property binSize Bin size in bps.
 * @property chromosomesSizes A map of chromosome name -> chromosome length entries.
 * @property dataQuery A query that returns a dataframe for each chromosome to serve as model input.
 * @property id A unique string identifier (include some kind of object hash if you compress identifiers). It's used
 * to generate the model file name if it's not provided. [reduceIds] is a recommended way to implement this property.
 */
interface SpanFitInformation {

    val build: String
    val binSize: Int
    val chromosomesSizes: LinkedHashMap<String, Int>
    val dataQuery: Query<Chromosome, DataFrame>
    val id: String

    fun genomeQuery(): GenomeQuery = GenomeQuery(Genome[build, chromosomesSizes], *chromosomesSizes.keys.toTypedArray())

    fun checkGenome(genome: Genome) {
        check(this.build == genome.build) {
            "Wrong genome build, expected: ${this.build}, got: ${genome.build}"
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
        checkGenome(chromosome.genome)
        checkChromosome(chromosome)
        return chromosome.range.slice(binSize).mapToInt { it.startOffset }.toArray()
    }

    /**
     * Returns boundaries of the squashed dataframe region corresponding to the given [chromosome].
     *
     * See also: [merge] and [split].
     */
    fun getChromosomesIndices(chromosome: Chromosome): Pair<Int, Int> {
        checkChromosome(chromosome)
        val offsetsMap = offsetsMap()
        val index = chromosomesSizes.keys.sorted().indexOf(chromosome.name)
        return offsetsMap[index] to offsetsMap[index + 1]
    }

    /**
     * Returns boundaries of the squashed dataframe region corresponding to the given [chromosome].
     *
     * See also: [merge] and [split].
     */
    fun getChromosomesIndices(chromosome: String): Pair<Int, Int> {
        check(chromosome in chromosomesSizes) {
            "Missing chromosome in ${chromosomesSizes.keys.toList()}: $chromosome"
        }
        val offsetsMap = offsetsMap()
        val index = chromosomesSizes.keys.sorted().indexOf(chromosome)
        return offsetsMap[index] to offsetsMap[index + 1]
    }

    /**
     * Merges (row-binds) the chromosome-wise dataframes in an unambiguous way.
     *
     * Inverse of [split].
     *
     * @param statesDataFrame a map of chromosome name -> dataframe entries. Must contain a dataframe with
     * row number equal to the number of bins on the appropriate chromosome for each chromosome in [chromosomesSizes].
     */
    fun merge(statesDataFrame: Map<String, DataFrame>): DataFrame {
        return DataFrame.rowBind(chromosomesSizes.keys.sorted().map { statesDataFrame[it]!! }.toTypedArray())
    }

    /**
     * Splits the squashed dataframe into chromosome-wise parts.
     *
     * Inverse of [merge].
     *
     * @param genomeQuery Optional smaller [GenomeQuery]. If provided, only requested chromosome-wise dataframes
     * are returned.
     */
    fun split(dataFrame: DataFrame, genomeQuery: GenomeQuery?): Map<String, DataFrame> {
        return if (genomeQuery != null) {
            checkGenome(genomeQuery.genome)
            genomeQuery.get()
                    .filter { it.name in chromosomesSizes }
                    .map { chromosome ->
                        val (start, end) = getChromosomesIndices(chromosome)
                        chromosome.name to dataFrame.iloc[start until end]
                    }.toMap()
        } else {
            chromosomesSizes.keys
                    .map { chromosome ->
                        val (start, end) = getChromosomesIndices(chromosome)
                        chromosome to dataFrame.iloc[start until end]
                    }.toMap()
        }
    }

    /**
     * Save the [SpanFitInformation] object at the given path as JSON.
     *
     * Inverse of [load].
     */
    fun save(path: Path) {
        path.parent.createDirectories()
        path.bufferedWriter().use { GSON.toJson(this, it) }
    }

    /**
     * Generates chromosome-wise dataframes for peak value calculation.
     *
     * If the map doesn't contain a specific chromosome, its peak values will be 0.0, so empty map is a perfectly
     * acceptable return value for this method.
     */
    fun scoresDataFrame(): Map<Chromosome, DataFrame>

    companion object {

        /**
         * Generate [SpanFitInformation.chromosomesSizes] instance from a [GenomeQuery]
         */
        fun chromSizes(genomeQuery: GenomeQuery) =
                LinkedHashMap<String, Int>().apply {
                    genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
                }

        object FragmentTypeAdapter : JsonSerializer<Fragment>, JsonDeserializer<Fragment> {

            override fun serialize(
                    src: Fragment, typeOfSrc: Type,
                    context: JsonSerializationContext
            ): JsonElement = context.serialize(src.toString())

            override fun deserialize(
                    json: JsonElement, typeOfT: Type,
                    context: JsonDeserializationContext
            ): Fragment {
                val str = context.deserialize<String>(json, object : TypeToken<String>() {}.type)
                try {
                    return Fragment.fromString(str)
                } catch (e: NumberFormatException) {
                    throw IllegalStateException("Failed to deserialize $str", e)
                }
            }
        }

        object PathTypeAdapter : JsonSerializer<Path>, JsonDeserializer<Path> {

            override fun serialize(
                    src: Path, typeOfSrc: Type,
                    context: JsonSerializationContext
            ): JsonElement = JsonPrimitive(src.toString())

            override fun deserialize(
                    json: JsonElement, typeOfT: Type,
                    context: JsonDeserializationContext
            ): Path {
                try {
                    return json.asString.toPath()
                } catch (e: NumberFormatException) {
                    throw IllegalStateException("Failed to deserialize ${json.asString}", e)
                }
            }
        }

        val GSON = GsonBuilder().setPrettyPrinting().setFieldNamingStrategy(
                    GSONUtil.NO_MY_UNDESCORE_NAMING_STRATEGY
                ).registerTypeAdapterFactory(
                    GSONUtil.classSpecificFactory(SpanFitInformation::class.java) { gson, factory ->
                        GSONUtil.classAndVersionAdapter(
                            gson, factory, "fit.information.fqn", "version"
                        )
                    }
                ).registerTypeAdapter(
                    object : TypeToken<Fragment>() {}.type, FragmentTypeAdapter
                ).registerTypeHierarchyAdapter(
                    Path::class.java, PathTypeAdapter
                ).create()

        /**
         * Loads a [SpanFitInformation] instance from a JSON file.
         *
         * Inverse of [SpanFitInformation.save]. Since "save" stores the fully-qualified class name,
         * "load" instantiates the correct class. If this class is not castable to [T],
         * [IllegalStateException] is thrown.
         */
        @Suppress("unchecked_cast")
        fun <T: SpanFitInformation> load(path: Path): T {
            val info = path.bufferedReader().use {
                GSON.fromJson(it, SpanFitInformation::class.java) as T?
            }
            check(info != null) { "Failed to load fit information from $path." }
            return info
        }
    }
}

/**
 * Contains the results of a Span-like model-fitting experiment.
 *
 * @property fitInfo The [SpanFitInformation] instance that describes the experiment input.
 * @property model The [ClassificationModel] that was fitted during the experiment.
 * @property logNullMemberships The chromosome-wise dataframes of log null probabilities, i.e.
 * the log probability of each observation under the null hypothesis. Each dataframe should at least contain
 * a column of floats or doubles labelled [SpanModelFitExperiment.NULL].
 */
data class SpanFitResults(
        val fitInfo: SpanFitInformation,
        val model: ClassificationModel,
        val logNullMemberships: Map<String, DataFrame>
) {
    companion object {
        internal val LOG = Logger.getLogger(SpanFitResults::class.java)
    }

    /**
     * @return Information about fit results including model and other parameters
     */
    fun about(): Map<String, String> {
        return when (model) {
            is MLFreeNBHMM -> {
                val signalMean = model.means[1]
                val noiseMean = model.means[0]
                mapOf(
                    "Signal mean" to signalMean.toString(),
                    "Noise mean" to noiseMean.toString(),
                    "Signal to noise" to ((signalMean + 1e-10) / (noiseMean + 1e-10)).toString()
                )
            }
            is PoissonRegressionMixture -> mapOf("Signal to noise" to model.signalToNoise.toString())
            else -> emptyMap()
        }
    }
}


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
        out Model : ClassificationModel, out FitInfo: SpanFitInformation, State : Any
> protected constructor(
        val fitInformation: FitInfo,
        private val modelFitter: Fitter<Model>,
        private val modelClass: Class<out Model>,
        private val availableStates: Array<State>,
        private val nullHypothesis: NullHypothesis<State>,
        private val fixedModelPath: Path? = null
) : Experiment("fit") {

    val genomeQuery = fitInformation.genomeQuery()
    val dataQuery = fitInformation.dataQuery

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
        return modelFitter.fit(preprocessedData, title = dataQuery.id)
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
                SpanFitInformation.load<SpanFitInformation>(informationPath)
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
                check(fitInfo == fitInformation) {
                    "Wrong model information: expected $fitInformation, but got $fitInfo"
                }
            }
        }
    }

    companion object {
        private const val INFORMATION_JSON = "information.json"
        private const val MODEL_JSON = "model.json"
        private const val NULL_NPZ = "null.npz"
        const val NULL = "null"

        val LOG: Logger = Logger.getLogger(ModelFitExperiment::class.java)

        /**
         * Retain only the chromosomes for which at least one treatment file has at least one read on them.
         */
        fun effectiveGenomeQuery(
                genomeQuery: GenomeQuery,
                paths: List<SpanDataPaths>,
                fragment: Fragment,
                unique: Boolean = true
        ): GenomeQuery {
            val chromosomes = genomeQuery.get()
            val nonEmptyChromosomes = hashSetOf<Chromosome>()
            paths.forEach { (t, c) ->
                val coverage = ReadsQuery(genomeQuery, t, unique, fragment, logFragmentSize = false).get()
                if (c != null) {
                    // we have to be sure that the control coverage cache is calculated for the full genome query,
                    // otherwise we can get some very hard-to-catch bugs later
                    ReadsQuery(genomeQuery, c, unique, fragment, logFragmentSize = false).get()
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
                // Sanity check
                genomeQuery?.let { info.checkGenome(it.genome) }
                val model = ClassificationModel.load<ClassificationModel>(dir / MODEL_JSON)
                val logNullMembershipsDF = DataFrame.load(dir / NULL_NPZ)
                val logNullMembershipsMap = info.split(logNullMembershipsDF, genomeQuery)

                LOG.info("Completed loading model: $tarPath")
                return@withTempDirectory SpanFitResults(info, model, logNullMembershipsMap)
            }
        }
    }
}

data class SpanDataPaths(
        val treatment: Path,
        val control: Path?
)

interface SpanAnalyzeFitInformation : SpanFitInformation {
    val data: List<SpanDataPaths>
    val fragment: Fragment
    val unique: Boolean
}
