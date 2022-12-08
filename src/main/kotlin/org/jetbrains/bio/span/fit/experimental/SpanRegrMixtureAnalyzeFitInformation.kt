package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.genome.sequence.CpGContent
import org.jetbrains.bio.span.coverage.NormalizedCoverageQuery
import org.jetbrains.bio.span.fit.AbstractSpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import org.jetbrains.bio.viktor.asF64Array
import java.nio.file.Path

data class SpanRegrMixtureAnalyzeFitInformation constructor(
    override val build: String,
    override val data: List<SpanDataPaths>,
    val mapabilityPath: Path?,
    override val fragment: Fragment,
    override val unique: Boolean,
    override val binSize: Int,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : AbstractSpanAnalyzeFitInformation {
    constructor(
        genomeQuery: GenomeQuery,
        data: SpanDataPaths,
        mapabilityPath: Path?,
        fragment: Fragment,
        unique: Boolean,
        binSize: Int
    ) : this(
        genomeQuery.build,
        listOf(data),
        mapabilityPath,
        fragment,
        unique,
        binSize,
        SpanFitInformation.chromSizes(genomeQuery)
    )

    override val id
        get() = reduceIds(
            listOfNotNull(data.single().treatment, data.single().control, mapabilityPath).map { it.stemGz } +
                    listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
        )

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            val genomeQuery = genomeQuery()
            val datum = data.single()

            return object : CachingQuery<Chromosome, DataFrame>() {

                private val treatmentCoverage = ReadsQuery(
                    genomeQuery, datum.treatment, unique, fragment, showLibraryInfo = false
                )
                private val controlCoverage = datum.control?.let {
                    ReadsQuery(genomeQuery, it, unique, fragment, showLibraryInfo = false)
                }

                override val id: String
                    get() = reduceIds(listOfNotNull(treatmentCoverage.id, controlCoverage?.id))

                override fun getUncached(input: Chromosome): DataFrame {
                    val y = binnedCoverage(input, treatmentCoverage.get(), binSize)
                    val control = controlCoverage?.let {
                        binnedCoverage(input, it.get(), binSize).map { it.toDouble() }.toDoubleArray()
                    }
                    val gc = CpGContent.binnedMeanCG(input, binSize)
                    val gc2 = (gc.asF64Array() * gc.asF64Array()).data
                    val mapability = mapabilityPath?.let { binnedMapability(input, it, binSize) }
                    var df = DataFrame().with("y", y)
                    df = df.with("GC", gc).with("GC2", gc2)
                    if (control != null) df = df.with("input", control)
                    if (mapability != null) df = df.with("mapability", mapability)
                    return df
                }
            }
        }

    @Transient
    var normalizedCoverageQuery: NormalizedCoverageQuery? = null

    @Synchronized
    override fun prepareData() {
        if (normalizedCoverageQuery == null) {
            normalizedCoverageQuery =
                NormalizedCoverageQuery(
                    genomeQuery(), data.single().treatment, data.single().control,
                    fragment, unique, binSize, showLibraryInfo = false
                )
        }
    }

    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(normalizedCoverageQuery != null) {
            "Please use prepareData before!"
        }
        return normalizedCoverageQuery!!.apply(chromosomeRange).toDouble()
    }

    override fun scaledTreatmentCoverage(chromosomeRange: ChromosomeRange): Double = 0.0
    override fun scaledControlCoverage(chromosomeRange: ChromosomeRange): Double? = null

    override fun hasControlData(): Boolean {
        return false
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION = 4

        fun binnedCoverage(chr: Chromosome, coverage: Coverage, binSize: Int): IntArray {
            return chr.range.slice(binSize).mapToInt { range ->
                coverage.getBothStrandsCoverage(range.on(chr))
            }.toArray()
        }

        fun binnedMapability(chr: Chromosome, mapabilityPath: Path, binSize: Int): DoubleArray {
            val bwFile = BigWigFile.read(mapabilityPath)
            // Mean genome mapability for all bins
            val meanMapability = bwFile.totalSummary.sum / bwFile.totalSummary.count
            return chr.range.slice(binSize).mapToDouble { range ->
                if (!bwFile.chromosomes.containsValue(chr.name)) {
                    return@mapToDouble meanMapability
                }
                bwFile.summarize(chr.name, range.startOffset, range.endOffset).single().sum / binSize
            }.toArray()
        }
    }
}