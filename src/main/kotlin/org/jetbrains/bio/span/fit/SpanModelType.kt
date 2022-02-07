package org.jetbrains.bio.span.fit

import org.jetbrains.bio.span.SpanCLA.LOG

enum class SpanModelType(val id: String, val extension: String, val description: String) {
    // Default SPAN model
    NB2Z_HMM("nb2zhmm", "span", "Negative binomial HMM 2 states with zero inflation"),
    NB2Z_MIXTURE("nb2zm", "span-nb2zm", "Negative binomial 2 states mixture with zero inflation"),
    NB3Z_HMM("nb3zhmm", "span-nb3zhmm", "Negative binomial HMM 3 states with zero inflation"),
    NB2_HMM("nb2hmm", "span-nb2hmm", "Negative binomial HMM 2 states"),
    NB3_HMM3("nb3hmm", "span-nb3hmm", "Negative binomial HMM 3 states"),
    POISSON_REGRESSION_MIXTURE(
        "p2zrm", "span-p2zrm", "Poisson 2 states regression mixture with zero inflation"
    ),
    NEGBIN_REGRESSION_MIXTURE(
        "nb2zrm", "span-nb2zrm", "Negative binomial 2 states regression mixture with zero inflation "
    );

    override fun toString() = description

    companion object {

        fun guessSpanModelById(id: String): SpanModelType {
            val model = values().find { it.id == id }
            check(model != null) {
                "Unrecognized value for --model-type command line option: $id"
            }
            return model
        }


        fun guessSpanModelByExtension(extension: String): SpanModelType? {
            val model = values().find { it.extension == extension }
            if (model == null) {
                LOG.warn("Unrecognized model extension '.$extension'")
            }
            return model
        }

    }
}
