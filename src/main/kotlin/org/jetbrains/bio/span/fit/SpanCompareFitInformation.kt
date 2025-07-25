package org.jetbrains.bio.span.fit

import com.google.common.math.DoubleMath
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.span.coverage.NormalizedCoverageQuery
import org.jetbrains.bio.span.coverage.binnedCoverageDataFrame
import org.jetbrains.bio.util.deleteIfExists
import org.jetbrains.bio.util.reduceIds

data class SpanCompareFitInformation(
    override val build: String,
    val data1: List<SpanDataPaths>,
    val data2: List<SpanDataPaths>,
    val labels1: List<String>,
    val labels2: List<String>,
    val explicitFormat: ReadsFormat?,
    val fragment: Fragment,
    val unique: Boolean,
    override val binSize: Int,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : SpanFitInformation {

    override val id
        get() = SpanAnalyzeFitInformation.generateId(data1 + data2, fragment, binSize, unique)

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            prepareData()
            val queries = normalizedCoverageQueries1!! + normalizedCoverageQueries2!!
            val labels = labels1 + labels2
            return object : CachingQuery<Chromosome, DataFrame>() {
                override fun getUncached(input: Chromosome): DataFrame {
                    return queries.binnedCoverageDataFrame(input, binSize, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(
                        queries.zip(labels).flatMap { (s, l) -> listOf(s.id, l) }
                                + listOf(binSize.toString())
                    )
            }
        }

    @Transient
    var normalizedCoverageQueries1: List<NormalizedCoverageQuery>? = null

    @Transient
    var normalizedCoverageQueries2: List<NormalizedCoverageQuery>? = null

    override fun prepareData() {
        if (normalizedCoverageQueries1 == null) {
            normalizedCoverageQueries1 = data1.map {
                NormalizedCoverageQuery(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    explicitFormat,
                    fragment,
                    unique,
                    binSize,
                    showLibraryInfo = false
                )
            }
        }
        if (normalizedCoverageQueries2 == null) {
            normalizedCoverageQueries2 = data2.map {
                NormalizedCoverageQuery(
                    genomeQuery(),
                    it.treatment,
                    it.control,
                    explicitFormat,
                    fragment,
                    unique,
                    binSize,
                    showLibraryInfo = false
                )
            }
        }
    }

    override fun isControlAvailable(): Boolean {
        return true
    }

    /**
     * Return log2 fold change of average summary coverage across data
     */
    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(normalizedCoverageQueries1 != null && normalizedCoverageQueries2 != null) {
            "Please use prepareData before!"
        }

        return if (normalizedCoverageQueries1!!.all { it.areCachesPresent() } &&
            normalizedCoverageQueries2!!.all { it.areCachesPresent() }) {
            val score1 = normalizedCoverageQueries1!!.sumOf { it.apply(chromosomeRange) }
                .toDouble() / normalizedCoverageQueries1!!.size
            val score2 = normalizedCoverageQueries2!!.sumOf { it.apply(chromosomeRange) }
                .toDouble() / normalizedCoverageQueries2!!.size
            if (score2 != 0.0) DoubleMath.log2(score1) - DoubleMath.log2(score2) else Double.MAX_VALUE
        } else {
            0.0
        }
    }

    override fun controlScore(chromosomeRange: ChromosomeRange): Double {
        check(normalizedCoverageQueries1 != null && normalizedCoverageQueries2 != null) {
            "Please use prepareData before!"
        }
        return normalizedCoverageQueries2!!.sumOf { it.score(chromosomeRange) } /
                normalizedCoverageQueries2!!.size
    }

    override fun cleanCaches() {
        normalizedCoverageQueries1?.forEach {
            it.treatmentReads.npzPath().deleteIfExists()
            it.controlReads?.npzPath()?.deleteIfExists()
        }
        normalizedCoverageQueries2?.forEach {
            it.treatmentReads.npzPath().deleteIfExists()
            it.controlReads?.npzPath()?.deleteIfExists()
        }
    }

    override fun toString(): String {
        return """
            SpanCompareFitInformation
                build=$build,
                data1=$data1,
                data2=$data2,
                fragment=$fragment,
                unique=$unique,
                binSize=$binSize,
                chromosomesSizes=$chromosomesSizes
        """.trimIndent()
    }

    override fun difference(loadedFitInfo: SpanFitInformation): String? {
        if (loadedFitInfo !is SpanCompareFitInformation) {
            return "Incompatible fit information type: ${loadedFitInfo::class.java.simpleName} " +
                    "instead of ${this::class.java.simpleName}"
        }
        if (binSize != loadedFitInfo.binSize) {
            return "Incompatible bin size: $binSize vs ${loadedFitInfo.binSize}"
        }
        if (fragment != loadedFitInfo.fragment) {
            return "Incompatible fragment: $fragment vs ${loadedFitInfo.fragment}"
        }
        if (genomeQuery().id != loadedFitInfo.genomeQuery().id) {
            return "Incompatible genome: ${genomeQuery().id} vs ${loadedFitInfo.genomeQuery().id}"
        }
        if (data1 != loadedFitInfo.data1) {
            return "Incompatible data1: $data1 vs ${loadedFitInfo.data1}"
        }
        if (data2 != loadedFitInfo.data2) {
            return "Incompatible data2: $data2 vs ${loadedFitInfo.data2}"
        }
        return null
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION: Int = 5

        fun effective(
            genomeQuery: GenomeQuery,
            paths1: List<SpanDataPaths>,
            paths2: List<SpanDataPaths>,
            labels1: List<String>,
            labels2: List<String>,
            explicitFormat: ReadsFormat?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
        ): SpanCompareFitInformation {
            return SpanCompareFitInformation(
                genomeQuery.build, paths1, paths2, labels1, labels2, explicitFormat, fragment, unique, binSize,
                SpanFitInformation.chromSizes(
                    SpanModelFitExperiment.filterGenomeQueryWithData(
                        genomeQuery, paths1 + paths2, explicitFormat, fragment, unique
                    )
                )
            )
        }
    }
}