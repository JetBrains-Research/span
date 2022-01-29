package org.jetbrains.bio.span.fit

enum class SpanModelType(val id: String, val extension: String, val description: String) {
    // Default SPAN model
    NB2Z_HMM("nb2zhmm", "span", "Negative binomial HMM 2states with zero inflation"),

    NB3Z_HMM("nb3zhmm", "span-nb3z", "Negative binomial HMM 3states with zero inflation"),
    NB2_HMM("nb2hmm", "span-nb2", "Negative binomial HMM 2states"),
    NB3_HMM3("nb3hmm", "span-nb3", "Negative binomial HMM 3states"),
    POISSON_REGRESSION_MIXTURE("p2zrm", "span-p2zrm", "Poisson regression mixture"),
    NEGBIN_REGRESSION_MIXTURE("nb2zrm", "span-nb2zrm", "Negative Binomial Regression mixture");

    override fun toString() = description

    companion object {

        fun guessSpanModelById(id: String): SpanModelType {
            val model = values().find { it.id == id }
            check(model != null) {
                "Unrecognized value for --model-type command line option: $id"
            }
            return model
        }


        fun guessSpanModelByExtension(extension: String): SpanModelType {
            val model = values().find { it.extension == extension }
            check(model != null) {
                "Unrecognized model extension '.$extension'"
            }
            return model
        }

    }
}
