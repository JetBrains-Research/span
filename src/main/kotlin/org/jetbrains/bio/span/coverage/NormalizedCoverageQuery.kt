package org.jetbrains.bio.span.coverage

import org.apache.commons.math3.stat.correlation.PearsonsCorrelation
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.span.fit.SpanPeakCallingExperiment.Companion.SPAN_DEFAULT_BIN
import org.jetbrains.bio.util.isAccessible
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import kotlin.math.abs
import kotlin.math.ceil
import kotlin.math.max

/**
 * Normalized coverage query.
 *
 * normCov = treatmentCov * treatmentScale - beta * controlCov * controlScale
 * treatmentScale and controlScale are computed to upscale smaller library.
 * beta is estimated as value between 0 and 1 minimizing absolute correlation between
 * normalized coverage and control coverage.
 * |correlation(treatmentCov * treatmentScale - beta * controlCov * controlScale, controlCov)|
 *
 * Scaling is inspired by https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137
 * Beta is inspired by https://www.cell.com/biophysj/fulltext/S0006-3495(17)30032-2#sec2

 * Reads are computed using [unique] option
 * Preprocessed reads are converted into tags using [Fragment] option
 */
class NormalizedCoverageQuery(
    val genomeQuery: GenomeQuery,
    val treatmentPath: Path,
    val controlPath: Path?,
    val fragment: Fragment,
    val unique: Boolean = true,
    val binSize: Int = SPAN_DEFAULT_BIN,
    val showLibraryInfo: Boolean = true,
) : Query<ChromosomeRange, Int> {

    override val id: String
        get() = reduceIds(
            listOfNotNull(
                treatmentPath.stemGz, controlPath?.stemGz, fragment.nullableInt, if (!unique) "keepdup" else null
            ).map { it.toString() }
        )

    override val description: String
        get() = "Treatment: $treatmentPath, Control: $controlPath, " +
                "Fragment: $fragment, Keep-dup: ${!unique}"

    val treatmentReads by lazy {
        ReadsQuery(genomeQuery, treatmentPath, unique, fragment, showLibraryInfo = showLibraryInfo)
    }

    val controlReads by lazy {
        controlPath?.let {
            ReadsQuery(genomeQuery, it, unique, fragment, showLibraryInfo = showLibraryInfo)
        }
    }

    /**
     * Shows whether the relevant caches are present.
     */
    val ready: Boolean
        get() = treatmentReads.npzPath().isAccessible() && controlReads?.npzPath()?.isAccessible() ?: true

    /**
     * Cached value for treatment and control scales, see [analyzeCoverage]
     */
    val coveragesNormalizedInfo by lazy {
        analyzeCoverage(genomeQuery, treatmentReads.get(), controlReads?.get(), binSize)
    }

    fun scaledTreatment(chromosomeRange: ChromosomeRange): Double {
        return treatmentReads.get().getBothStrandsCoverage(chromosomeRange) * coveragesNormalizedInfo.treatmentScale
    }

    fun scaledControl(chromosomeRange: ChromosomeRange): Double? {
        if (controlReads == null) {
            return null
        }
        return controlReads!!.get().getBothStrandsCoverage(chromosomeRange) * coveragesNormalizedInfo.controlScale
    }

    override fun apply(t: ChromosomeRange): Int {
        if (controlPath == null) {
            return treatmentReads.get().getBothStrandsCoverage(t)
        }
        val (treatmentScale, controlScale, beta) = coveragesNormalizedInfo
        val treatmentCoverage = treatmentReads.get().getBothStrandsCoverage(t)
        val controlCoverage = controlReads!!.get().getBothStrandsCoverage(t)
        return max(
            0,
            ceil(treatmentCoverage * controlScale - controlCoverage * treatmentScale * beta).toInt()
        )
    }

    companion object {
        private val LOG: Logger = LoggerFactory.getLogger(NormalizedCoverageQuery::class.java)

        data class NormalizedCoverageInfo(
            val treatmentScale: Double,
            val controlScale: Double,
            val beta: Double
        )


        /**
         * Compute coefficients to normalize coverage, see [NormalizedCoverageInfo].
         * Whenever possible use [coveragesNormalizedInfo] cached value
         */
        fun analyzeCoverage(
            genomeQuery: GenomeQuery,
            treatmentCoverage: Coverage,
            controlCoverage: Coverage?,
            binSize: Int
        ): NormalizedCoverageInfo {
            if (controlCoverage == null) {
                return NormalizedCoverageInfo(1.0, 1.0, 0.0)
            }
            val treatmentTotal = genomeQuery.get().sumOf {
                treatmentCoverage.getBothStrandsCoverage(it.range.on(it)).toLong()
            }
            val controlTotal = genomeQuery.get().sumOf {
                controlCoverage.getBothStrandsCoverage(it.range.on(it)).toLong()
            }
            val rationTreatmentToControl = 1.0 * treatmentTotal / controlTotal
            val controlScale: Double
            val treatmentScale: Double
            // Scale to the biggest depth
            if (rationTreatmentToControl <= 1) {
                treatmentScale = 1.0 / rationTreatmentToControl
                controlScale = 1.0
                LOG.info(
                    "Upscale treatment($treatmentTotal) to control($controlTotal) x${"%.3f".format(treatmentScale)}"
                )
            } else {
                controlScale = rationTreatmentToControl
                treatmentScale = 1.0
                LOG.info(
                    "Upscale control($controlTotal) to treatment($treatmentTotal) x${"%.3f".format(controlScale)}"
                )
            }
            LOG.debug("Estimating beta")
            val beta = estimateBeta(
                genomeQuery, treatmentCoverage, treatmentScale, controlCoverage, controlScale, binSize
            )
            LOG.info("Beta for signal control correction: ${"%.3f".format(beta)}")
            return NormalizedCoverageInfo(treatmentScale, controlScale, beta)
        }

        private fun estimateBeta(
            genomeQuery: GenomeQuery,
            treatmentCoverage: Coverage,
            treatmentScale: Double,
            controlCoverage: Coverage,
            controlScale: Double,
            bin: Int = SPAN_DEFAULT_BIN,
            betaStep: Double = 0.01
        ): Double {
            val maxChromosome = genomeQuery.get().maxByOrNull { it.length }!!
            val binnedTreatment = maxChromosome.range.slice(bin).mapToDouble { range ->
                treatmentCoverage.getBothStrandsCoverage(range.on(maxChromosome)) * treatmentScale
            }.toArray()
            val binnedControl = maxChromosome.range.slice(bin).mapToDouble { range ->
                controlCoverage.getBothStrandsCoverage(range.on(maxChromosome)) * controlScale
            }.toArray()
            var b = 0.0
            var minCorrelation = 1.0
            var minB = -1.0
            // Reuse array to reduce GC
            val binnedNorm = DoubleArray(binnedTreatment.size)
            val pearsonCorrelation = PearsonsCorrelation()
            while (b < 1) {
                for (i in binnedNorm.indices) {
                    binnedNorm[i] = binnedTreatment[i] - b * binnedControl[i]
                }
                val c = abs(pearsonCorrelation.correlation(binnedNorm, binnedControl))
                if (c < minCorrelation) {
                    minCorrelation = c
                    minB = b
                }
                b += betaStep
            }
            return minB
        }
    }

}

fun List<NormalizedCoverageQuery>.binnedCoverageDataFrame(
    chromosome: Chromosome,
    binSize: Int,
    labels: Array<String>
): DataFrame {
    var res = DataFrame()
    forEachIndexed { d, inputQuery ->
        val binnedCoverage = chromosome.range.slice(binSize).mapToInt { range ->
            inputQuery.apply(range.on(chromosome))
        }.toArray()
        res = res.with(labels[d], binnedCoverage)
    }
    return res
}