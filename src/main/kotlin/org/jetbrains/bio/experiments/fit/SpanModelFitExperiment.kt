package org.jetbrains.bio.experiments.fit

import com.google.common.annotations.VisibleForTesting
import com.google.common.math.IntMath
import com.google.gson.*
import com.google.gson.reflect.TypeToken
import kotlinx.support.jdk7.use
import org.apache.log4j.Logger
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanModelFitExperiment.Companion.createEffectiveQueries
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.sequence.TwoBitSequence
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
import java.lang.reflect.Type
import java.math.RoundingMode
import java.nio.file.Path
import java.util.*
import kotlin.collections.LinkedHashMap

/**
 * Since all the chromosomes are squashed in [SpanModelFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 *
 * [labels] refer to the coverage dataframe column labels, not to the supervised learning annotations.
 */
data class SpanFitInformation(
        val build: String,
        val data: List<SpanPathsToData>,
        val labels: List<String>,
        val fragment: Fragment,
        val unique: Boolean,
        val binSize: Int,
        val chromosomesSizes: LinkedHashMap<String, Int>,
        val version: Int
) {

    constructor(
            genomeQuery: GenomeQuery,
            paths: List<SpanPathsToData>,
            labels: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
    ): this(
        genomeQuery.build,
        paths,
        labels, fragment, unique, binSize,
        LinkedHashMap<String, Int>().apply {
            genomeQuery.get().sortedBy { it.name }.forEach { this[it.name] = it.length }
        },
        VERSION
    )

    fun genomeQuery(): GenomeQuery = GenomeQuery(Genome[build, chromosomesSizes], *chromosomesSizes.keys.toTypedArray())

    internal fun checkGenome(genome: Genome) {
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

    fun scoresDataFrame(): Map<Chromosome, DataFrame> {
        val gq = genomeQuery()
        val queries = data.map {
            CoverageScoresQuery(gq, it.pathTreatment, it.pathInput, fragment, binSize, unique)
        }
        if (queries.any { !it.ready }) {
            return emptyMap()
        }
        return gq.get().associateBy({it}) {
            queries.scoresDataFrame(it, labels.toTypedArray())
        }
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

    internal fun getChromosomesIndices(chromosome: String): Pair<Int, Int> {
        check(chromosome in chromosomesSizes) {
            "Missing chromosome in ${chromosomesSizes.keys.toList()}: $chromosome"
        }
        val offsetsMap = offsetsMap()
        val index = chromosomesSizes.keys.sorted().indexOf(chromosome)
        return offsetsMap[index] to offsetsMap[index + 1]
    }

    internal fun save(path: Path) {
        path.parent.createDirectories()
        path.bufferedWriter().use { GSON.toJson(this, it) }
    }

    internal fun merge(statesDataFrame: Map<String, DataFrame>): DataFrame {
        return DataFrame.rowBind(chromosomesSizes.keys.sorted().map { statesDataFrame[it]!! }.toTypedArray())
    }

    internal fun split(dataFrame: DataFrame, genomeQuery: GenomeQuery?): Map<String, DataFrame> {
        if (genomeQuery != null) {
            checkGenome(genomeQuery.genome)
            return genomeQuery.get()
                    .filter { it.name in chromosomesSizes }
                    .map { chromosome ->
                        val (start, end) = getChromosomesIndices(chromosome)
                        chromosome.name to dataFrame.iloc[start until end]
                    }.toMap()
        } else {
            return chromosomesSizes.keys
                    .map { chromosome ->
                        val (start, end) = getChromosomesIndices(chromosome)
                        chromosome to dataFrame.iloc[start until end]
                    }.toMap()
        }
    }

    companion object {
        const val VERSION: Int = 2

        /**
         * Using Treatment and Control class instead of [Pair] here for nice GSON serialization
         */
        data class TC(val path: String, val control: String?)

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
            ): JsonElement = context.serialize(src.toString())

            override fun deserialize(
                    json: JsonElement, typeOfT: Type,
                    context: JsonDeserializationContext
            ): Path {
                val str = context.deserialize<String>(json, object : TypeToken<String>() {}.type)
                try {
                    return str.toPath()
                } catch (e: NumberFormatException) {
                    throw IllegalStateException("Failed to deserialize $str", e)
                }
            }
        }


        private val GSON = GsonBuilder()
                .registerTypeAdapter(object : TypeToken<Fragment>() {}.type, FragmentTypeAdapter)
                .registerTypeAdapter(object : TypeToken<Path>() {}.type, PathTypeAdapter)
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
            else -> emptyMap()
        }
    }
}


/**
 * A generic class for Span (Semi-supervised Peak Analyzer) - tool for analyzing and comparing ChIP-Seq data.
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
        paths: List<SpanPathsToData>,
        labels: List<String>,
        fragment: Fragment,
        val binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        availableStates: Array<State>,
        private val nullHypothesis: NullHypothesis<State>,
        val unique: Boolean = true,
        private val fixedModelPath: Path? = null
) : ModelFitExperiment<Model, State>(
    createEffectiveQueries(externalGenomeQuery, paths, labels, fragment, binSize, unique),
    modelFitter, modelClass, availableStates
) {
    fun getIntCover(chr1: Chromosome, coverage: Coverage, bin: Int): IntArray {
        val len = (chr1.length - 1) / bin + 1
        val cover = IntArray(len)
        for (i in 0 until len - 1) {
            cover[i] = coverage.getBothStrandsCoverage(ChromosomeRange(i * bin, (i + 1) * bin, chr1))
        }
        cover[len - 1] = coverage.getBothStrandsCoverage(ChromosomeRange((len-1) * bin, chr1.length, chr1))
        return cover
    }
    fun getDoubleCover(chr1: Chromosome, coverage: Coverage, bin: Int): DoubleArray {
        val len = (chr1.length - 1) / bin + 1
        val cover = DoubleArray(len)
        for (i in 0 until len - 1) {
            cover[i] = coverage
                    .getBothStrandsCoverage(ChromosomeRange(i * bin, (i + 1) * bin, chr1))
                    .toDouble()
        }
        cover[len - 1] = coverage
                .getBothStrandsCoverage(ChromosomeRange((len-1) * bin, chr1.length, chr1))
                .toDouble()
        return cover
    }
    fun getGC(chr1: Chromosome, bin: Int): DoubleArray {
        val len = (chr1.length - 1) / bin + 1
        val seq: TwoBitSequence = chr1.sequence
        val GCcontent = DoubleArray(len)
        for (i in 0 until len - 1) {
            GCcontent[i] = seq.substring(i*bin, (i + 1)*bin).count { it == 'c' || it == 'g' }.toDouble()/bin
        }
        GCcontent[len - 1] = seq
                .substring((len-1)*bin, seq.length)
                .count { it == 'c'|| it == 'g' }
                .toDouble()/( seq.length - (len-1)*bin)
        return GCcontent
    }

    fun getMappability(chr1: Chromosome, path_mappability: Path, bin: Int): DoubleArray {
        if (BigWigFile.read(path_mappability).chromosomes.containsValue(chr1.name)) {
            val mapSummary = BigWigFile
                    .read(path_mappability)
                    .summarize(chr1.name, 0, chr1.length - chr1.length % bin, numBins = (chr1.length - 1) / bin)
            val result = DoubleArray(mapSummary.size + 1) {
                if (it < mapSummary.size) mapSummary[it].sum / bin else 1.0
            }
            result[mapSummary.size] = BigWigFile
                    .read(path_mappability)
                    .summarize(chr1.name, chr1.length - chr1.length % bin, 0)[0].sum / chr1.length % bin
            return result
        }
        val meanMappability = BigWigFile.read(path_mappability).totalSummary.sum/BigWigFile.read(path_mappability).totalSummary.count
        return DoubleArray(chr1.length) {meanMappability}
    }

    private val preprocessedData: List<Preprocessed<DataFrame>> by lazy {
        preprocessData(paths)
    }

    private fun preprocessData(paths: List<SpanPathsToData>): List<Preprocessed<DataFrame>> {
        val readsQueryTreatment = ReadsQuery(genomeQuery, paths[0].pathTreatment, unique)
        val readsQueryInput = paths[0].pathInput?.let { ReadsQuery(genomeQuery, it, unique) }
        val coverageTreatment = readsQueryTreatment.get()
        val coverageInput = readsQueryInput?.get()
        val chrList: List<Chromosome> = genomeQuery.get().sortedBy { it.name }

        val coverLength = chrList.map { it.length/binSize + 1 }.sum()
        val coverTreatment = IntArray (coverLength)
        val coverInput = DoubleArray (coverLength)
        val GCcontent = DoubleArray (coverLength)
        val mappability = DoubleArray (coverLength)
        var prevIdx = 0
        chrList.forEach {
            val arraySize = it.length / binSize + 1
            System.arraycopy(
                    getIntCover(it, coverageTreatment, binSize),
                    0, coverTreatment, prevIdx, arraySize)
            if (coverageInput != null) {
                System.arraycopy(
                    getDoubleCover(it, coverageInput, binSize),
                    0, coverInput, prevIdx, arraySize)
            }
            System.arraycopy(
                    getGC(it, binSize),
                    0, GCcontent, prevIdx, arraySize)
            if (paths[0].pathMappability != null) {
                System.arraycopy(
                        getMappability(it, paths[0].pathMappability!!, binSize),
                        0, mappability, prevIdx, arraySize)
            }
            prevIdx += (arraySize)
        }

        var covar = DataFrame().with("y", coverTreatment)
                .with("GC", GCcontent)
                .with("GC2", DoubleArray(GCcontent.size) { GCcontent[it] * GCcontent[it] })

        if (paths[0].pathMappability != null) {
            covar = covar.with("mappability", mappability)
        }

        if (coverageInput != null) {
            covar = covar.with("input", coverInput)
        }

        return listOf(Preprocessed.of(covar))
    }

    val results: SpanFitResults by lazy {
        getOrLoadResults()
    }

    override fun getStatesDataFrame(chromosome: Chromosome): DataFrame = sliceStatesDataFrame(statesDataFrame, chromosome)

    override fun doCalculations() {
        results.logNullMemberships
    }

    val fitInformation = SpanFitInformation(genomeQuery, paths, labels, fragment, unique, binSize)




    // XXX It is important to use get() here, because id is overridden in superclasses
    private val modelPath: Path
        get() = fixedModelPath ?: experimentPath / "$id.span"

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


        @VisibleForTesting
        /**
         * Create pair of
         * 1. Effective genomeQuery, i.e. only chromosomes with some reads on them
         * 2. Data query required for [ModelFitExperiment]
         */
        internal fun createEffectiveQueries(
                genomeQuery: GenomeQuery,
                paths: List<SpanPathsToData>,
                labels: List<String>,
                fragment: Fragment,
                binSize: Int,
                unique: Boolean = true
        ): Pair<GenomeQuery, Query<Chromosome, DataFrame>> {
            val chromosomes = genomeQuery.get()
            val nonEmptyChromosomes = hashSetOf<Chromosome>()
            paths.forEach { (t, _) ->
                val coverage = ReadsQuery(
                    genomeQuery, t,
                    unique = unique, fragment = fragment, logFragmentSize = false
                ).get()
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
            val effectiveGenomeQuery = GenomeQuery(
                genomeQuery.genome,
                *nonEmptyChromosomes.map { it.name }.toTypedArray()
            )
            return effectiveGenomeQuery to object : CachingQuery<Chromosome, DataFrame>() {
                val scores = paths.map {
                    CoverageScoresQuery(genomeQuery, it.pathTreatment, it.pathInput, fragment, binSize, unique)
                }

                override fun getUncached(input: Chromosome): DataFrame {
                    return scores.scoresDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(scores.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
            }
        }


        fun loadResults(genomeQuery: GenomeQuery? = null, tarPath: Path): SpanFitResults {
            LOG.info("Loading model: $tarPath")
            return withTempDirectory(tarPath.stem) { dir ->
                LOG.debug("Started model file decompress: $tarPath")
                Tar.decompress(tarPath, dir.toFile())

                LOG.debug("Completed model file decompress and started loading: $tarPath")
                val info = SpanFitInformation.load(dir / INFORMATION_JSON)
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

data class SpanPathsToData(
        val pathTreatment: Path,
        val pathInput: Path?,
        val pathMappability: Path?
) {
    constructor(pathTreatment: Path, pathInput: Path)
            : this(pathTreatment, pathInput, null)
}