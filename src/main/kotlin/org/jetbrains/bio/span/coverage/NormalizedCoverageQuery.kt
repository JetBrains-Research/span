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
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BETA_STEP
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_BIN
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
 * A class representing a normalized coverage query.
 *
 * normCov = treatmentCov - controlCov * controlScale * beta
 *
 * beta is estimated as value between 0 and 1 minimizing absolute correlation between
 *      normalized coverage and control coverage:
 *      | correlation(treatmentCov - controlCov * controlScale * beta, controlCov) |
 *
 * Scaling is inspired by https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137
 * Beta is inspired by https://www.cell.com/biophysj/fulltext/S0006-3495(17)30032-2#sec2
 *
 * @property genomeQuery The genome query.
 * @property treatmentPath The path to the treatment file.
 * @property controlPath The path to the control file.
 * @property fragment The fragment.
 * @property unique Flag indicating whether to use unique reads.
 * @property binSize The bin size.
 * @property showLibraryInfo Flag indicating whether to show library information.
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
                "Fragment: $fragment, Keep-duplicates: ${!unique}"

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
    fun areCachesPresent(): Boolean =
        treatmentReads.npzPath().isAccessible() && (controlReads?.npzPath()?.isAccessible() ?: true)

    /**
     * Cached value for treatment and control scales, see [analyzeCoverage]
     */
    val coveragesNormalizedInfo by lazy {
        analyzeCoverage(genomeQuery, treatmentReads.get(), controlReads?.get(), binSize)
    }

    fun score(t: ChromosomeRange): Double {
        return treatmentReads.get().getBothStrandsCoverage(t).toDouble()
    }

    fun isControlAvailable(): Boolean {
        return controlReads != null
    }

    fun controlScore(chromosomeRange: ChromosomeRange): Double {
        require (controlReads != null) { "Control is not available" }
        return controlReads!!.get().getBothStrandsCoverage(chromosomeRange) * coveragesNormalizedInfo.controlScale
    }

    override fun apply(t: ChromosomeRange): Int {
        val treatmentCoverage = treatmentReads.get().getBothStrandsCoverage(t)
        if (controlPath == null) {
            return treatmentCoverage
        }
        val (controlScale, beta, _) = coveragesNormalizedInfo
        val controlCoverage = controlReads!!.get().getBothStrandsCoverage(t)
        return max(0.0, ceil(treatmentCoverage - controlCoverage * controlScale * beta)).toInt()
    }

    companion object {
        private val LOG: Logger = LoggerFactory.getLogger(NormalizedCoverageQuery::class.java)

        data class NormalizedCoverageInfo(
            val controlScale: Double,
            val beta: Double,
            val minCorrelation: Double
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
                return NormalizedCoverageInfo(0.0, 0.0, 0.0)
            }
            val treatmentTotal = genomeQuery.get().sumOf {
                treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange).toLong()
            }
            val controlTotal = genomeQuery.get().sumOf {
                controlCoverage.getBothStrandsCoverage(it.chromosomeRange).toLong()
            }
            // Scale control to treatment
            val controlScale = 1.0 * treatmentTotal  / controlTotal
            val (beta, minCorrelation) = estimateBeta(
                genomeQuery, treatmentCoverage, controlCoverage, controlScale, binSize
            )
            LOG.info(
                "Treatment ${"%,d".format(treatmentTotal)}, " +
                        "control ${"%,d".format(controlTotal)} x ${"%.3f".format(controlScale)}, " +
                        "min correlation ${"%.3f".format(minCorrelation)}, beta ${"%.3f".format(beta)}"
            )
            return NormalizedCoverageInfo(controlScale, beta, minCorrelation)
        }

        private fun estimateBeta(
            genomeQuery: GenomeQuery,
            treatmentCoverage: Coverage,
            controlCoverage: Coverage,
            controlScale: Double,
            bin: Int,
            betaStep: Double = SPAN_DEFAULT_BETA_STEP,
        ): Pair<Double, Double> {
            // Estimate beta corrected signal only on not empty chromosomes
            val chromosomeWithMaxSignal = genomeQuery.get()
                .maxByOrNull { treatmentCoverage.getBothStrandsCoverage(it.chromosomeRange) } ?:
                return 0.0 to 0.0
            val binnedScaledTreatment = chromosomeWithMaxSignal.range.slice(bin).mapToDouble {
                treatmentCoverage.getBothStrandsCoverage(it.on(chromosomeWithMaxSignal)).toDouble()
            }.toArray()
            val binnedScaledControl = chromosomeWithMaxSignal.range.slice(bin).mapToDouble {
                controlCoverage.getBothStrandsCoverage(it.on(chromosomeWithMaxSignal)) * controlScale
            }.toArray()
            var b = 0.0
            var minCorrelation = 1.0
            var minB = 0.0
            // Reuse array to reduce GC
            val binnedNorm = DoubleArray(binnedScaledTreatment.size)
            val pearsonCorrelation = PearsonsCorrelation()
            while (b < 1) {
                for (i in binnedNorm.indices) {
                    binnedNorm[i] = binnedScaledTreatment[i] - binnedScaledControl[i] * b
                }
                val c = abs(pearsonCorrelation.correlation(binnedNorm, binnedScaledControl))
                if (c <= minCorrelation) {
                    minCorrelation = c
                    minB = b
                }
                b += betaStep
            }
            if (minB == 0.0) {
                LOG.warn("Failed to estimate beta-value for control correction")
            }
            return minB to minCorrelation
        }
    }

}

/**
 * Calculates the binned coverage DataFrame for a list of normalized coverage queries.
 * Using clip percentile may slightly help to avoid out-of-range during model fit.
 */
fun List<NormalizedCoverageQuery>.binnedCoverageDataFrame(
    chromosome: Chromosome,
    binSize: Int,
    labels: Array<String>,
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