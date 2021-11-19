package org.jetbrains.bio.experiments.fit

import com.google.common.math.IntMath
import com.google.gson.*
import com.google.gson.reflect.TypeToken
import kotlinx.support.jdk7.use
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.statistics.gson.GSONUtil
import org.jetbrains.bio.util.*
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
                    gson, factory,
                    "fit.information.fqn", "version"
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
        fun <T : SpanFitInformation> load(path: Path): T {
            val info = path.bufferedReader().use {
                GSON.fromJson(it, SpanFitInformation::class.java) as T?
            }
            check(info != null) { "Failed to load fit information from $path." }
            return info
        }
    }
}