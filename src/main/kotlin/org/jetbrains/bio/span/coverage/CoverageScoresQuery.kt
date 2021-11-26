package org.jetbrains.bio.span.coverage

import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.coverage.Coverage
import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.query.Query
import org.jetbrains.bio.genome.query.ReadsQuery
import org.jetbrains.bio.util.exists
import org.jetbrains.bio.util.reduceIds
import org.jetbrains.bio.util.stemGz
import java.nio.file.Path
import kotlin.math.ceil
import kotlin.math.max
import kotlin.math.min

/**
 * Coverage query to reproduce DiffBind-like scores.
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

    private val treatmentCoverage = ReadsQuery(
        genomeQuery, treatmentPath, unique, fragment, showLibraryInfo = showLibraryInfo
    )

    private val controlCoverage = controlPath?.let {
        ReadsQuery(genomeQuery, it, unique, fragment, showLibraryInfo = showLibraryInfo)
    }

    /**
     * Shows whether the relevant caches are present.
     */
    val ready: Boolean
        get() = treatmentCoverage.npzPath().exists && controlCoverage?.npzPath()?.exists ?: true

    val scale: Double = if (controlCoverage != null)
        computeScale(genomeQuery, treatmentCoverage.get(), controlCoverage.get())
    else
        0.0

    override fun apply(t: ChromosomeRange): Int {
        return getScore(t, treatmentCoverage.get(), controlCoverage?.get(), scale)
    }

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

        internal fun getScore(
            t: ChromosomeRange,
            treatmentReads: Coverage,
            controlReads: Coverage?,
            scale: Double
        ): Int {
            return if (controlReads != null) {
                val plusStrandBin = t.on(Strand.PLUS)
                val plusStrandScore = max(
                    0,
                    treatmentReads.getCoverage(plusStrandBin) -
                            ceil(controlReads.getCoverage(plusStrandBin) * scale).toInt()
                )
                val minusStrandBin = t.on(Strand.MINUS)
                val minusStrandScore = max(
                    0,
                    treatmentReads.getCoverage(minusStrandBin) -
                            ceil(controlReads.getCoverage(minusStrandBin) * scale).toInt()
                )
                plusStrandScore + minusStrandScore
            } else {
                treatmentReads.getBothStrandsCoverage(t)
            }
        }
    }

}