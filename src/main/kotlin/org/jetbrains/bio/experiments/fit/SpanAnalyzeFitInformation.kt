package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.genome.coverage.Fragment

interface SpanAnalyzeFitInformation : SpanFitInformation {
    val data: List<SpanDataPaths>
    val fragment: Fragment
    val unique: Boolean
}
