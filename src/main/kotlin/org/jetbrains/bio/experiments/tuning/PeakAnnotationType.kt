package org.jetbrains.bio.experiments.tuning

import java.awt.Color

enum class PeakAnnotationType(val id: String,
                              val desc: String,
                              val help: String,
                              val color: Color) {
    // Use colors from original markup (see article)
    PEAK_START("peakStart",
               "Peak Start","Range contains only one peak start",
               Color(111, 188, 112)),
    PEAK_END("peakEnd",
             "Peak End","Range contains only one peak end",
             Color(255, 76, 76)),
    NO_PEAKS("noPeaks",
             "No Peaks", "Range doesn't contain peaks",
             Color.lightGray),
    PEAKS("peaks",
          "Peaks","Range contains at least 1 peak",
          Color(164, 69, 238));

    override fun toString() = id

    companion object {
        fun from(name: String) = values().firstOrNull { it.id == name }
    }
}