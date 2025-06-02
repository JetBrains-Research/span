package org.jetbrains.bio.span.fit

import org.jetbrains.bio.genome.coverage.Fragment
import org.jetbrains.bio.genome.format.ReadsFormat

interface AbstractSpanAnalyzeFitInformation : SpanFitInformation {
    val paths: List<SpanDataPaths>
    val explicitFormat: ReadsFormat?
    val fragment: Fragment
    val unique: Boolean
}
