package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.span.coverage.NormalizedCoverageQuery
import org.jetbrains.bio.span.coverage.binnedCoverageDataFrame
import org.jetbrains.bio.util.deleteIfExists
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz

/**
 * Since all the chromosomes are squashed in [SpanModelFitExperiment] and processed by the single model,
 * this class is used to access chromosomes information from that model.
 *
 * See [getChromosomesIndices] and [offsets] for details.
 *
 * [labels] refer to the coverage dataframe column labels, not to the supervised learning annotations.
 */
data class SpanAnalyzeFitInformation(
    override val build: String,
    override val paths: List<SpanDataPaths>,
    val labels: List<String>,
    override val fragment: Fragment,
    override val unique: Boolean,
    override val binSize: Int,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : AbstractSpanAnalyzeFitInformation {

    constructor(
        genomeQuery: GenomeQuery,
        paths: List<SpanDataPaths>,
        labels: List<String>,
        fragment: Fragment,
        unique: Boolean,
        binSize: Int
    ) : this(
        genomeQuery.build, paths,
        labels, fragment, unique, binSize,
        SpanFitInformation.chromSizes(genomeQuery)
    )

    override val id: String
        get() = generateId(paths, fragment, binSize, unique)


    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            prepareData()
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return normalizedCoverageQueries!!.binnedCoverageDataFrame(
                        input, binSize, labels.toTypedArray()
                    )
                }

                override val id: String
                    get() = reduceIds(
                        normalizedCoverageQueries!!.zip(labels)
                            .flatMap { (s, l) -> listOf(s.id, l) } + listOf(binSize.toString())
                    )
            }
        }

    @Transient
    var normalizedCoverageQueries: List<NormalizedCoverageQuery>? = null

    override fun prepareData() {
        if (normalizedCoverageQueries == null) {
            normalizedCoverageQueries = paths.map {
                NormalizedCoverageQuery(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    fragment,
                    unique,
                    binSize,
                    showLibraryInfo = true
                ).apply {
                    // Force computing caches
                    treatmentReads.get()
                    controlReads?.get()
                }
            }
        }
    }

    /**
     * Returns average coverage by tracks
     */
    override fun score(chromosomeRange: ChromosomeRange): Double {
        return normalizedCoverageQueries!!.sumOf { it.score(chromosomeRange) } /
                normalizedCoverageQueries!!.size
    }

    override fun isControlAvailable(): Boolean =
        normalizedCoverageQueries!!.all { it.controlReads != null }

    override fun controlScore(chromosomeRange: ChromosomeRange): Double {
        check(isControlAvailable()) {
            "Control is not available"
        }
        return normalizedCoverageQueries!!.sumOf { it.controlScore(chromosomeRange) } /
                normalizedCoverageQueries!!.size
    }

    override fun cleanCaches() {
        normalizedCoverageQueries?.forEach {
            it.treatmentReads.npzPath().deleteIfExists()
            it.controlReads?.npzPath()?.deleteIfExists()
        }
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 5

        /**
         * Generates a model ID based on the provided parameters.
         *
         * @param paths A list of [SpanDataPaths] representing treatment and control pairs.
         * @param binSize The bin size.
         * @param fragment The fragment type.
         * @param unique Indicates whether the experiment is unique.
         * @return The generated ID
         */
        fun generateId(
            paths: List<SpanDataPaths>,
            fragment: Fragment,
            binSize: Int,
            unique: Boolean
        ) = reduceIds(
            paths.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOfNotNull(fragment.nullableInt, binSize).map { it.toString() } +
                    listOfNotNull(if (unique) "unique" else null)
        )

        fun createFitInformation(
            genomeQuery: GenomeQuery,
            paths: List<SpanDataPaths>,
            labels: List<String>,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
        ): SpanAnalyzeFitInformation {
            val genomeQueryWithData =
                SpanModelFitExperiment.filterGenomeQueryWithData(genomeQuery, paths, fragment, unique)
            return SpanAnalyzeFitInformation(
                genomeQueryWithData.build,
                paths,
                labels,
                fragment,
                unique,
                binSize,
                SpanFitInformation.chromSizes(genomeQueryWithData)
            )
        }

        /**
         * Optimization to avoid synchronized lazy on NormalizedCoverageQuery#treatmentReads
         * Replacing calls NormalizedCoverageQuery#score with non-blocked direct realization
         */
        fun fastScore(treatmentCoverages: List<Coverage>, range: ChromosomeRange) =
            treatmentCoverages.sumOf { it.getBothStrandsCoverage(range) }.toDouble() / treatmentCoverages.size

        /**
         * Optimization to avoid synchronized lazy on NormalizedCoverageQuery#controlReads
         * Replacing calls NormalizedCoverageQuery#controlScore with non-blocked direct realization
         */
        fun fastControlScore(
            controlCoverages: List<Coverage>,
            controlScales: List<Double>,
            chromosomeRange: ChromosomeRange
        ) = controlCoverages.zip(controlScales)
            .sumOf { (c, s) -> c.getBothStrandsCoverage(chromosomeRange) * s } / controlCoverages.size


    }
}
