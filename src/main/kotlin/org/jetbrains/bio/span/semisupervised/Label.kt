package org.jetbrains.bio.span.semisupervised

enum class Label(
    val id: String,
    val desc: String,
    val help: String,
) {

    PEAKS("peaks", "Peaks", "Range contains at least 1 peak"),
    NO_PEAKS("noPeaks", "No Peaks", "Range doesn't contain peaks"),
    PEAK_START("peakStart", "Peak Start", "Range contains only one peak start"),
    PEAK_END("peakEnd", "Peak End", "Range contains only one peak end");

    override fun toString() = id

    companion object {
        fun from(name: String) = values().firstOrNull { it.id == name }
    }
}