package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.genome.containers.GenomeMap

data class SpanPeaksResult(
    val fdr: Double,
    val sensitivity: Double,
    val sensitivitySummits: Double?,
    val gap: Int,
    val peaks: GenomeMap<List<Peak>>
) {
    fun toList(): List<Peak> = peaks.genomeQuery.get().flatMap { peaks[it] }
}