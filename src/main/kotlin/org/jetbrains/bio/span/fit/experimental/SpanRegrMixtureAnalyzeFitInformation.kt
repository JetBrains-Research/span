package org.jetbrains.bio.span.fit.experimental

import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.genome.sequence.CpGContent
import org.jetbrains.bio.span.coverage.NormalizedCoverageQuery
import org.jetbrains.bio.span.fit.AbstractSpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation.Companion.generateId
import org.jetbrains.bio.span.fit.SpanDataPaths
import org.jetbrains.bio.span.fit.SpanFitInformation
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.viktor.asF64Array
import java.nio.file.Path

data class SpanRegrMixtureAnalyzeFitInformation(
    override val build: String,
    override val paths: List<SpanDataPaths>,
    override val explicitFormat: ReadsFormat?,
    val mapabilityPath: Path?,
    override val fragment: Fragment,
    override val unique: Boolean,
    override val binSize: Int,
    override val chromosomesSizes: LinkedHashMap<String, Int>
) : AbstractSpanAnalyzeFitInformation {
    constructor(
        genomeQuery: GenomeQuery,
        data: SpanDataPaths,
        explicitFormat: ReadsFormat?,
        mapabilityPath: Path?,
        fragment: Fragment,
        unique: Boolean,
        binSize: Int
    ) : this(
        genomeQuery.build,
        listOf(data),
        explicitFormat,
        mapabilityPath,
        fragment,
        unique,
        binSize,
        SpanFitInformation.chromSizes(genomeQuery)
    )

    override val id
        get() = generateId(paths, fragment, binSize, unique)

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            val genomeQuery = genomeQuery()
            val datum = paths.single()

            return object : CachingQuery<Chromosome, DataFrame>() {

                private val treatmentCoverage = ReadsQuery(
                    genomeQuery, datum.treatment, explicitFormat, unique, fragment, showLibraryInfo = false
                )
                private val controlCoverage = datum.control?.let {
                    ReadsQuery(genomeQuery, it, explicitFormat, unique, fragment, showLibraryInfo = false)
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

    override fun prepareData() {
        if (normalizedCoverageQuery == null) {
            normalizedCoverageQuery =
                NormalizedCoverageQuery(
                    genomeQuery(), paths.single().treatment, paths.single().control,
                    explicitFormat, fragment, unique, binSize, showLibraryInfo = true
                )
        }
    }

    override fun isControlAvailable(): Boolean {
        return normalizedCoverageQuery?.controlReads != null &&
                normalizedCoverageQuery!!.controlReads!!.isAccessible()
    }

    override fun score(chromosomeRange: ChromosomeRange): Double {
        check(normalizedCoverageQuery != null) {
            "Please use prepareData before!"
        }
        return normalizedCoverageQuery!!.apply(chromosomeRange).toDouble()
    }

    override fun controlScore(chromosomeRange: ChromosomeRange): Double = TODO("Not implemented")

    override fun cleanCaches() {
        // Not implemented yet
    }

    override fun difference(loadedFitInfo: SpanFitInformation): String? {
        // Not yet implemented
        return null
    }

    companion object {
        @Suppress("MayBeConstant", "unused")
        @Transient
        @JvmField
        val VERSION = 5

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