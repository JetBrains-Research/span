package org.jetbrains.bio.span.fit

import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.statistics.model.ClassificationModel
import org.slf4j.LoggerFactory

/**
 * Contains the extended results of a Span-like model-fitting experiment.
 * See [SpanFitResults]
 */
class SpanFitResultsExt(
    fitInfo: SpanFitInformation,
    model: ClassificationModel,
    logNullMemberships: Map<String, DataFrame>,
    val coveragesDataFrameMap: Map<String, DataFrame>,
    val statesDataFrameMap: Map<String, DataFrame>,
) : SpanFitResults(fitInfo, model, logNullMemberships) {
    companion object {
        internal val LOG = LoggerFactory.getLogger(SpanFitResultsExt::class.java)
    }
}
