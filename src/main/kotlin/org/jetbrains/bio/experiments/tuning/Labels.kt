package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.big.ExtendedBedEntry
import org.jetbrains.bio.dataset.CellId
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.dataset.ReplicateDataKey
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.genome.Location
import org.jetbrains.bio.genome.LocationAware
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.io.BedField
import org.jetbrains.bio.io.BedFormat
import org.jetbrains.bio.io.toLocation
import org.jetbrains.bio.io.unpackRegularFields
import org.jetbrains.bio.util.toPath
import java.nio.file.Path

data class LabelledTrack(val cellId: CellId,
                         val name: String,
                         val trackPath: Path,
                         val labelPath: Path,
                         val failed: Boolean,
                         val testLabelPath: Path? = null)


private val LABELS_KEY = ReplicateDataKey("labels") {
    it?.toString()?.toPath()
}

private val TEST_LABELS_KEY = ReplicateDataKey("test_labels") {
    it?.toString()?.toPath()
}

/**
 * Uses custom label path provided in the "labels" metainfo field of the replicate data.\
 * In absence of this, uses the label path provided in the "labels" field of auxiliary data.
 * In absence of both, the track is discarded.
 */
fun DataConfig.extractLabelledTracks(target: String): List<LabelledTrack> {
    val auxLabelPath = this[target, LABELS_KEY]
    return tracksMap.entries
            .filter { it.key.dataType == target }
            .flatMap { entry -> entry.value.map { it to entry.key } }
            .filter { auxLabelPath != null || it.first.second[LABELS_KEY] != null }
            .map { (data, key) ->
                val customLabelPath = data.second[LABELS_KEY]
                val testLabelPath = data.second[TEST_LABELS_KEY] ?: this[target, TEST_LABELS_KEY]
                LabelledTrack(key.cellId, data.first, data.second.path,
                        customLabelPath ?: auxLabelPath!!, data.second.failedTrack, testLabelPath)
            }
}


data class PeakAnnotation(override val location: Location,
                          val type: PeakAnnotationType) : LocationAware, Comparable<PeakAnnotation> {

    override fun compareTo(other: PeakAnnotation): Int = compareValuesBy(this, other, { it.location }, { it.type })

    /**
     * Checks according to [peaks] provided
     */
    fun check(peaks: LocationsMergingList): Boolean {
        val intersection = peaks.intersect(location)
        return when (type) {
            PeakAnnotationType.NO_PEAKS -> intersection.isEmpty()
            PeakAnnotationType.PEAKS -> intersection.isNotEmpty()
            PeakAnnotationType.PEAK_START -> intersection.singleOrNull()
                    ?.let { it.startOffset != location.startOffset } == true
            PeakAnnotationType.PEAK_END -> intersection.singleOrNull()
                    ?.let { it.endOffset != location.endOffset } == true
        }
    }

    fun asBedEntry() = ExtendedBedEntry(location.chromosome.name,
                                        location.startOffset, location.endOffset,
                                        type.id, itemRgb = type.color.rgb)

    companion object {
        fun loadLabels(labelPath: Path, genome: Genome): List<PeakAnnotation> {
            val format = BedFormat(BedField.NAME)
            return format.parse(labelPath) { it.toList() }.map {
                val e = it.unpackRegularFields(format)
                val type = PeakAnnotationType.from(e.name)
                requireNotNull("Unknown type: $type")
                PeakAnnotation(e.toLocation(genome), type!!)
            }
        }
    }
}

