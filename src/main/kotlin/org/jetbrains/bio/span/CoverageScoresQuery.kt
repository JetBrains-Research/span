package org.jetbrains.bio.span

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.GenomeMap
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.CachingQuery
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.util.exists
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import java.nio.file.Path
import kotlin.math.ceil
import kotlin.math.max
import kotlin.math.min

/**
 * Coverage query to reproduce DiffBind-like scores for binned genome.
 * 1. Reads are computed using [unique] option
 * 2. Preprocessed reads are converted into tags using [Fragment] option
 * 3. Control coverage is scaled to match treatment
 * 4. Resulting score is treatment - scale * control
 */
class CoverageScoresQuery(
    val genomeQuery: GenomeQuery,
    private val treatmentPath: Path,
    private val controlPath: Path?,
    val fragment: Fragment,
    val binSize: Int,
    val unique: Boolean = true
) : CachingQuery<Chromosome, IntArray>() {


    override val id: String
        get() = reduceIds(
            listOfNotNull(
                treatmentPath.stemGz, controlPath?.stemGz, fragment.nullableInt,
                binSize, if (!unique) "keepdup" else null
            ).map { it.toString() }
        )

    override val description: String
        get() = "Treatment: $treatmentPath, Control: $controlPath, " +
                "Fragment: $fragment, Bin: $binSize, Keep-dup: ${!unique}"

    override fun getUncached(input: Chromosome): IntArray {
        return scores[input]
    }

    private val treatmentCoverage = ReadsQuery(
        genomeQuery, treatmentPath, unique = unique, fragment = fragment
    )

    private val controlCoverage = controlPath?.let {
        ReadsQuery(genomeQuery, it, unique = unique, fragment = fragment)
    }

    val scores: GenomeMap<IntArray> by lazy {
        val treatmentCoverage = treatmentCoverage.get()
        val controlCoverage = controlCoverage?.get()
        val scale: Double = if (controlCoverage != null)
            computeScale(genomeQuery, treatmentCoverage, controlCoverage)
        else
            0.0
        return@lazy genomeMap(genomeQuery, parallel = true) {
            computeBinnedScores(it, treatmentCoverage, controlCoverage, binSize, scale)
        }
    }

    /**
     * Shows whether the relevant caches are present.
     */
    val ready: Boolean
        get() = treatmentCoverage.npzPath().exists && controlCoverage?.npzPath()?.exists ?: true

    companion object {

        /**
         * Compute multiplier to scale control coverage to match treatment
         */
        internal fun computeScale(
            genomeQuery: GenomeQuery,
            treatmentCoverage: Coverage,
            controlCoverage: Coverage
        ): Double {
            val conditionSize = genomeQuery.get().sumOf {
                treatmentCoverage.getCoverage(it.range.on(it).on(Strand.PLUS)).toLong() +
                        treatmentCoverage.getCoverage(it.range.on(it).on(Strand.MINUS)).toLong()
            }
            val controlSize = genomeQuery.get().sumOf {
                controlCoverage.getCoverage(it.range.on(it).on(Strand.PLUS)).toLong() +
                        controlCoverage.getCoverage(it.range.on(it).on(Strand.MINUS)).toLong()
            }
            return min(1.0, conditionSize.toDouble() / controlSize)
        }

        /**
         * Returns the scores of a given [chromosome] sliced into binSizes with [binSize] width.
         * Score is precisely [Coverage] within bins if no control given, or DiffBind-like score.
         */
        internal fun computeBinnedScores(
            chromosome: Chromosome,
            treatmentCoverage: Coverage,
            controlCoverage: Coverage?,
            binSize: Int,
            scale: Double
        ): IntArray {
            return chromosome.range.slice(binSize).mapToInt { bin ->
                val plusStrandBin = bin.on(chromosome, Strand.PLUS)
                val minusStrandBin = bin.on(chromosome, Strand.MINUS)
                if (controlCoverage != null) {
                    val plusStrandScore = max(
                        0,
                        treatmentCoverage.getCoverage(plusStrandBin) -
                                ceil(controlCoverage.getCoverage(plusStrandBin) * scale).toInt()
                    )
                    val minusStrandScore = max(
                        0,
                        treatmentCoverage.getCoverage(minusStrandBin) -
                                ceil(controlCoverage.getCoverage(minusStrandBin) * scale).toInt()
                    )
                    return@mapToInt plusStrandScore + minusStrandScore
                } else {
                    return@mapToInt treatmentCoverage.getBothStrandsCoverage(bin.on(chromosome))
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
