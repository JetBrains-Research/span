package org.jetbrains.bio.span.semisupervised

import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.LocationAware
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.format.BedField
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.toLocation
import org.jetbrains.bio.genome.format.unpackRegularFields
import java.nio.file.Path

data class LocationLabel(
    override val location: Location,
    val type: Label
) : LocationAware, Comparable<LocationLabel> {

    override fun compareTo(other: LocationLabel): Int = compareValuesBy(this, other, { it.location }, { it.type })

    /**
     * Checks according to [peaks] provided
     */
    fun check(peaks: LocationsMergingList): Boolean {
        val intersection = peaks.intersect(location)
        return when (type) {
            Label.NO_PEAKS -> intersection.isEmpty()
            Label.PEAKS -> intersection.isNotEmpty()
            Label.PEAK_START -> intersection.singleOrNull()
                ?.let { it.startOffset != location.startOffset } == true

            Label.PEAK_END -> intersection.singleOrNull()
                ?.let { it.endOffset != location.endOffset } == true
        }
    }


    companion object {

        /**
         * Loads labels from path
         */
        fun loadLabels(labelPath: Path, genome: Genome): List<LocationLabel> {
            val format = BedFormat(BedField.NAME)
            return format.parse(labelPath) { it.toList() }.map {
                val e = it.unpackRegularFields(format)
                val type = Label.from(e.name)
                requireNotNull(type) { "Unknown type: ${e.name}" }
                LocationLabel(e.toLocation(genome), type)
            }
        }

        /**
         * Function to compute Peaks [peaks] vs Labels [labels]
         * @return map like structure
         */
        fun computeErrors(labels: List<LocationLabel>, peaks: LocationsMergingList): LabelErrors =
            LabelErrors().apply {
                labels.forEach { observe(it, it.check(peaks)) }
            }

    }
}

