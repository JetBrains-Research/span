package org.jetbrains.bio.span.coverage

import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.util.isAccessible
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.nio.file.Path
import kotlin.math.max

/**
 * Coverage query with scaling down to smallest of treatment or control libraries
 * 1. Reads are computed using [unique] option
 * 2. Preprocessed reads are converted into tags using [Fragment] option
 * 3. Treatment and control scales are computed
 */
class CoverageScoresQuery(
    val genomeQuery: GenomeQuery,
    private val treatmentPath: Path,
    private val controlPath: Path?,
    val fragment: Fragment,
    val unique: Boolean = true,
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
     * Cached value for treatment and control scales, see [computeScales]
     */
    val treatmentAndControlScales by lazy {
        computeScales(genomeQuery, treatmentReads.get(), controlReads?.get())
    }

    fun scaledTreatment(chromosomeRange: ChromosomeRange): Double {
        return treatmentReads.get().getBothStrandsCoverage(chromosomeRange) * treatmentAndControlScales.first
    }

    fun scaledControl(chromosomeRange: ChromosomeRange): Double? {
        if (controlReads == null) {
            return null
        }
        return controlReads!!.get().getBothStrandsCoverage(chromosomeRange) * treatmentAndControlScales.second
    }

    override fun apply(t: ChromosomeRange): Int {
        if (controlPath == null) {
            return treatmentReads.get().getBothStrandsCoverage(t)
        }
        val (treatmentScale, controlScale) = treatmentAndControlScales
        return max(
            0, (controlScale * treatmentReads.get().getBothStrandsCoverage(t) -
                    treatmentScale * controlReads!!.get().getBothStrandsCoverage(t)).toInt()
        )
    }

    companion object {
        private val LOG: Logger = LoggerFactory.getLogger(CoverageScoresQuery::class.java)

        /**
         * Compute multipliers to upscale smaller library to bigger,
         * whenever possible use [treatmentAndControlScales] cached value
         */
        fun computeScales(
            genomeQuery: GenomeQuery,
            treatmentCoverage: Coverage,
            controlCoverage: Coverage?
        ): Pair<Double, Double> {
            if (controlCoverage == null) {
                return 1.0 to 1.0
            }
            val treatmentTotal = genomeQuery.get().sumOf {
                treatmentCoverage.getBothStrandsCoverage(it.range.on(it)).toLong()
            }
            val controlTotal = genomeQuery.get().sumOf {
                controlCoverage.getBothStrandsCoverage(it.range.on(it)).toLong()
            }
            val rationTreatmentToControl = 1.0 * treatmentTotal / controlTotal
            val controlScale: Double
            val treatScale: Double
            // Scale to the biggest depth
            if (rationTreatmentToControl <= 1) {
                treatScale = 1.0 / rationTreatmentToControl
                controlScale = 1.0
                LOG.debug(
                    "Upscale treatment($treatmentTotal) to control($controlTotal) x${"%.3f".format(treatScale)}"
                )
            } else {
                controlScale = rationTreatmentToControl
                treatScale = 1.0
                LOG.debug(
                    "Upscale control($controlTotal) to treatment($treatmentTotal) x${"%.3f".format(controlScale)}"
                )
            }
            return treatScale to controlScale
        }
    }

}