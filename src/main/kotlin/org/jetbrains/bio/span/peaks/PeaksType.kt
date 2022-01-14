package org.jetbrains.bio.span.peaks

import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.util.CancellableState

enum class PeaksType(val cmd: String) {
    PEAKS_TYPE_ISLANDS("islands"),
    PEAKS_TYPE_FDR_GAP("simple");

    companion object {
        fun getTypeFromCmd(cmd: String): PeaksType = when (cmd) {
            PEAKS_TYPE_ISLANDS.cmd -> PEAKS_TYPE_ISLANDS
            PEAKS_TYPE_FDR_GAP.cmd -> PEAKS_TYPE_FDR_GAP
            else -> throw  IllegalStateException("Unknown type: $cmd")
        }

    }
}

fun SpanFitResults.getPeaks(
    genomeQuery: GenomeQuery,
    fdr: Double,
    gap: Int,
    cancellableState: CancellableState? = null,
    peaksType: PeaksType = PeaksType.PEAKS_TYPE_ISLANDS
): List<Peak> =
    when (peaksType) {
        PeaksType.PEAKS_TYPE_ISLANDS -> getIslands(genomeQuery, fdr, gap, cancellableState = cancellableState)
        PeaksType.PEAKS_TYPE_FDR_GAP -> getFdrGapPeaks(genomeQuery, fdr, gap, cancellableState = cancellableState)
    }

fun SpanFitResults.getPeaks(
    chromosome: Chromosome,
    fdr: Double,
    gap: Int,
    peaksType: PeaksType = PeaksType.PEAKS_TYPE_ISLANDS
): List<Peak> =
    when (peaksType) {
        PeaksType.PEAKS_TYPE_ISLANDS -> getChromosomeIslands(chromosome, fdr, gap)
        PeaksType.PEAKS_TYPE_FDR_GAP -> getChromosomeFdrGapPeaks(chromosome, fdr, gap)
    }
