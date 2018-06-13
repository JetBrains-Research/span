package org.jetbrains.bio.experiments.tuning

import org.jetbrains.bio.datasets.CellId
import org.jetbrains.bio.experiments.DataConfig
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.tools.Washu
import java.nio.file.Path

object MACS2 : Tool2Tune<Double>() {
    override val id = "macs_narrow"
    override val suffix = ".narrowPeak"
    const val DEFAULT_FDR = 0.05
    const val DEFAULT_ULI_FDR = 1E-4
    val FDRS = arrayListOf(1E-2, DEFAULT_FDR, DEFAULT_ULI_FDR, 1E-6, 1E-8, 1E-10)
    override val parameters = FDRS
    override val transform = Double::toString

    override fun callPeaks(configuration: DataConfig, p: Path, parameter: Double) {
        Washu().runMACS2(Genome[configuration.genome], p, parameter, false)
    }

    override fun fileName(cellId: CellId, replicate: String, target: String, parameter: Double): String {
        return "${cellId}_${replicate}_${target}_fdr${parameter}_peaks.narrowPeak"
    }

    override fun defaultParams(uli: Boolean) = if (uli) DEFAULT_ULI_FDR else DEFAULT_FDR
}

object MACS2_BROAD : Tool2Tune<Double>() {
    override val id = "macs_broad"
    override val suffix = ".broadPeak"

    const val DEFAULT_FDR = 0.1
    const val DEFAULT_ULI_FDR = 1E-4
    val FDRS = arrayListOf(DEFAULT_FDR, 1E-2, DEFAULT_ULI_FDR, 1E-6, 1E-8, 1E-10)
    override val parameters = FDRS
    override val transform = Double::toString

    override fun callPeaks(configuration: DataConfig, p: Path, parameter: Double) {
        Washu().runMACS2(Genome[configuration.genome], p, parameter, true)
    }

    override fun fileName(cellId: CellId, replicate: String, target: String, parameter: Double): String {
        return "${cellId}_${replicate}_${target}_fdr${parameter}_peaks.broadPeak"
    }

    override fun defaultParams(uli: Boolean) = if (uli) DEFAULT_ULI_FDR else DEFAULT_FDR
}


object SICER : Tool2Tune<Pair<Double, Int>>() {
    override val id = "sicer"
    override val suffix = "-island.bed"

    const val DEFAULT_FDR = 1E-2
    const val DEFAULT_ULI_FDR = 1E-6
    val FDRS = arrayListOf(DEFAULT_FDR, 1E-4, DEFAULT_ULI_FDR, 1E-8)
    const val DEFAULT_GAP = 600
    val GAPS = arrayListOf(0, 200, DEFAULT_GAP, 1200)
    override val parameters = FDRS.flatMap { fdr -> GAPS.map { gap -> fdr to gap } }
    override val transform: (Pair<Double, Int>) -> String = { (fdr, gap) -> "${fdr}_$gap" }

    override fun callPeaks(configuration: DataConfig, p: Path, parameter: Pair<Double, Int>) {
        Washu().runSicer(Genome[configuration.genome], p, parameter.first, listOf("200", "150", parameter.second.toString()))
    }

    override fun fileName(cellId: CellId, replicate: String, target: String, parameter: Pair<Double, Int>): String {
        return "${cellId}_${replicate}_$target-W200-G${parameter.second}-FDR${parameter.first}-island.bed"
    }

    override fun defaultParams(uli: Boolean) = if (uli)
        DEFAULT_ULI_FDR to DEFAULT_GAP else DEFAULT_FDR to DEFAULT_GAP
}