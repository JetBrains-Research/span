package org.jetbrains.bio.span.fit

import java.nio.file.Path

/**
 * Data class to represent treatment and control pair.
 * Later it is used to produce DiffBind like coverage, see [CoverageScoresQuery]
 */
data class SpanDataPaths(val treatment: Path, val control: Path?)

