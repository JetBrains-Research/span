package org.jetbrains.bio.experiments.fit.experimental

import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.experiments.fit.SpanDataPaths
import org.jetbrains.bio.experiments.fit.SpanFitInformation
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.genome.query.stemGz
import org.jetbrains.bio.genome.sequence.CpGContent
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.viktor.asF64Array
import java.nio.file.Path

data class Span2FitInformation constructor(
        override val build: String,
        override val data: List<SpanDataPaths>,
        val mapabilityPath: Path?,
        override val fragment: Fragment,
        override val unique: Boolean,
        override val binSize: Int,
        override val chromosomesSizes: LinkedHashMap<String, Int>
) : SpanAnalyzeFitInformation {
    constructor(
            genomeQuery: GenomeQuery,
            data: SpanDataPaths,
            mapabilityPath: Path?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
    ) : this(
            genomeQuery.build, listOf(data), mapabilityPath, fragment, unique, binSize, SpanFitInformation.chromSizes(genomeQuery)
    )

    override val id
        get() = reduceIds(
                listOfNotNull(data.single().treatment, data.single().control, mapabilityPath).map { it.stemGz } +
                        listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
        )

    override fun scoresDataFrame(): Map<Chromosome, DataFrame> {
        val datum = data.single()
        return genomeQuery().get().associateBy({ it }) {
            DataFrame().with("coverage", binnedCoverage(
                    it,
                    ReadsQuery(genomeQuery(), datum.treatment, unique, fragment, logFragmentSize = false).get(),
                    binSize
            ))
        }

    }

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            val genomeQuery = genomeQuery()
            val datum = data.single()

            return object : CachingQuery<Chromosome, DataFrame>() {

                private val treatmentCoverage = ReadsQuery(
                        genomeQuery, datum.treatment, unique, fragment, logFragmentSize = false
                )
                private val controlCoverage = datum.control?.let {
                    ReadsQuery(genomeQuery, it, unique, fragment, logFragmentSize = false)
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

    companion object {

        const val VERSION = 3

        fun binnedCoverage(chr: Chromosome, coverage: Coverage, binSize: Int): IntArray {
            val len = (chr.length - 1) / binSize + 1
            val res = IntArray(len)
            for (i in 0 until len - 1) {
                res[i] = coverage.getBothStrandsCoverage(
                        ChromosomeRange(i * binSize, (i + 1) * binSize, chr)
                )
            }
            res[len - 1] = coverage.getBothStrandsCoverage(
                    ChromosomeRange((len - 1) * binSize, chr.length, chr)
            )
            return res
        }

        fun binnedMapability(chr: Chromosome, mapabilityPath: Path, binSize: Int): DoubleArray {
            val bwFile = BigWigFile.read(mapabilityPath)

            val len = (chr.length - 1) / binSize + 1

            if (!bwFile.chromosomes.containsValue(chr.name)) {
                // the chromosome isn't present in the bigWig file, use mean genome mapability for all bins
                val meanMapability = bwFile.totalSummary.sum / bwFile.totalSummary.count
                return DoubleArray(len) { meanMapability }
            }

            val mapSummary = bwFile.summarize(chr.name, 0, (len - 1) * binSize, numBins = len - 1)
            val res = DoubleArray(len)
            for (i in 0 until len - 1) {
                res[i] = mapSummary[i].sum / binSize
            }
            res[len - 1] = bwFile.summarize(chr.name, (len - 1) * binSize, chr.length).single().sum /
                    (chr.length - (len - 1) * binSize)
            return res
        }
    }
}