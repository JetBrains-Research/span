package org.jetbrains.bio.span

import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.data.DataFrame
import java.nio.file.Path

class CoverageScoresQuery(val genomeQuery: GenomeQuery,
                          private val treatmentPath: Path,
                          private val controlPath: Path?,
                          val fragment: Int?,
                          val binSize: Int) : CachingQuery<Chromosome, IntArray>() {


    override val id: String
        get() = reduceIds(listOfNotNull(treatmentPath.stemGz, controlPath?.stemGz, fragment, binSize).map {
            it.toString()
        })

    override val description: String
        get() = "Treatment: $treatmentPath, Control: $controlPath, Fragment: $fragment, Bin: $binSize"

    override fun getUncached(input: Chromosome): IntArray {
        return scores[input]
    }

    val scores: GenomeMap<IntArray> by lazy {
        val treatmentCoverage = ReadsQuery(genomeQuery, treatmentPath, unique = true, fragment = fragment).get()
        val controlCoverage = if (controlPath != null)
            ReadsQuery(genomeQuery, controlPath, unique = true, fragment = fragment).get()
        else
            null
        val scale: Double = if (controlCoverage != null)
            computeScale(genomeQuery, treatmentCoverage, controlCoverage)
        else
            0.0
        return@lazy genomeMap(genomeQuery, parallel = true) {
            computeScores(it, treatmentCoverage, controlCoverage, binSize, scale)
        }
    }

    companion object {

        internal fun computeScale(genomeQuery: GenomeQuery,
                                  treatmentCoverage: Coverage,
                                  controlCoverage: Coverage): Double {
            val conditionSize = genomeQuery.get().map {
                treatmentCoverage.getCoverage(it.range.on(it).on(Strand.PLUS)).toLong() +
                        treatmentCoverage.getCoverage(it.range.on(it).on(Strand.MINUS)).toLong()
            }.sum()
            val controlSize = genomeQuery.get().map {
                controlCoverage.getCoverage(it.range.on(it).on(Strand.PLUS)).toLong() +
                        controlCoverage.getCoverage(it.range.on(it).on(Strand.MINUS)).toLong()
            }.sum()
            return Math.min(1.0, conditionSize.toDouble() / controlSize)
        }

        /**
         * Returns the scores of a given [chromosome] sliced into binSizes with [binSize] width.
         * Score is precisely [Coverage] within bins if no control given, or DiffBind-like score.
         */
        internal fun computeScores(chromosome: Chromosome,
                                   treatmentCoverage: Coverage,
                                   controlCoverage: Coverage?,
                                   binSize: Int,
                                   scale: Double): IntArray {
            return chromosome.range.slice(binSize).mapToInt { b ->
                val plusStrandBin = b.on(chromosome, Strand.PLUS)
                val minusStrandBin = b.on(chromosome, Strand.MINUS)
                if (controlCoverage != null) {
                    val plusStrandScore = Math.max(0,
                            treatmentCoverage.getCoverage(plusStrandBin) -
                                    Math.ceil(controlCoverage.getCoverage(plusStrandBin) * scale).toInt())
                    val minusStrandScore = Math.max(0,
                            treatmentCoverage.getCoverage(minusStrandBin) -
                                    Math.ceil(controlCoverage.getCoverage(minusStrandBin) * scale).toInt())
                    return@mapToInt plusStrandScore + minusStrandScore
                } else {
                    return@mapToInt treatmentCoverage.getBothStrandsCoverage(b.on(chromosome))
                }
            }.toArray()
        }

    }
}

fun List<CoverageScoresQuery>.scoresDataFrame(chromosome: Chromosome, labels: Array<String>): DataFrame {
    var res = DataFrame()
    forEachIndexed { d, inputQuery ->
        val binnedCoverage = inputQuery.apply(chromosome)
        res = res.with(labels[d], binnedCoverage)
    }
    return res
}
