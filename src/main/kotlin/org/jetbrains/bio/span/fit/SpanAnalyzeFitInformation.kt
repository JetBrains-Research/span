package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.span.coverage.BinnedNormalizedCoverageQuery
import org.jetbrains.bio.span.coverage.NormalizedCoverageQuery
import org.jetbrains.bio.span.coverage.binnedCoverageDataFrame
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
    override val data: List<SpanDataPaths>,
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
        get() = reduceIds(
            data.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                    listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
        )

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            prepareData()
            val binnedNormalizedCoverageQueries =
                normalizedCoverageQueries!!.map { BinnedNormalizedCoverageQuery(it, binSize) }
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return binnedNormalizedCoverageQueries.binnedCoverageDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(binnedNormalizedCoverageQueries.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
            }
        }

    @Transient
    var normalizedCoverageQueries: List<NormalizedCoverageQuery>? = null

    @Synchronized
    override fun prepareData() {
        if (normalizedCoverageQueries == null) {
            normalizedCoverageQueries = data.map {
                NormalizedCoverageQuery(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    fragment,
                    unique,
                    showLibraryInfo = false
                )
            }
        }
    }

    /**
     * Returns average coverage by tracks
     */
    @Synchronized
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        return if (normalizedCoverageQueries!!.all { it.ready }) {
            normalizedCoverageQueries!!.sumOf { it.apply(chromosomeRange) }
                .toDouble() / normalizedCoverageQueries!!.size
        } else {
            0.0
        }
    }

    @Synchronized
    override fun scaledTreatmentCoverage(chromosomeRange: ChromosomeRange): Double {
        check(normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        return normalizedCoverageQueries!!.sumOf { it.scaledTreatment(chromosomeRange) } /
                normalizedCoverageQueries!!.size
    }

    @Synchronized
    override fun scaledControlCoverage(chromosomeRange: ChromosomeRange): Double? {
        check(normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        if (normalizedCoverageQueries!!.any { it.controlReads == null }) {
            return null
        }
        return normalizedCoverageQueries!!.sumOf { it.scaledControl(chromosomeRange)!! } /
                normalizedCoverageQueries!!.size
    }

    override fun hasControlData(): Boolean {
        check(normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        return normalizedCoverageQueries!!.any { it.controlPath != null }
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 4

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
    }
}
