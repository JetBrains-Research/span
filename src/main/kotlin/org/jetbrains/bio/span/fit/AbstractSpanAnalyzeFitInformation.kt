package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.coverage.Fragment

interface AbstractSpanAnalyzeFitInformation : SpanFitInformation {
    val paths: List<SpanDataPaths>
    val fragment: Fragment
    val unique: Boolean
}
