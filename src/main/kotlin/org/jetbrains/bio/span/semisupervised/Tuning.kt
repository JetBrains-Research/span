package org.jetbrains.bio.span.semisupervised

import org.jetbrains.bio.dataframe.DataFrame
import org.slf4j.LoggerFactory
import java.nio.file.Path

class TuningResultTable {

    val names: MutableList<String> = arrayListOf()
    private val parameters: MutableList<String> = arrayListOf()
    private val errors: Array<MutableList<Double>> = Array(Label.values().size + 1) {
        arrayListOf()
    }

    fun minError() = errors.mapNotNull { it.minOrNull() }.minOrNull()

    fun addRecord(name: String, parameter: String, labelErrors: LabelErrors) {
        names.add(name)
        parameters.add(parameter)
        Label.values().forEachIndexed { i, type -> errors[i].add(labelErrors.error(type)) }
        errors[Label.values().size].add(labelErrors.error())
    }

    fun save(path: Path) {
        var dataFrame = DataFrame()
            .with("name", names.toTypedArray())
            .with("parameter", parameters.toTypedArray())
            .with("error", errors[Label.values().size].toDoubleArray())
        Label.values().forEachIndexed { i, type ->
            dataFrame = dataFrame.with("error_$type", errors[i].toDoubleArray())
        }
        dataFrame.save(path)
    }
}

class TuningResults(
    private val results: TuningResultTable = TuningResultTable(),
    private val optimalResults: TuningResultTable = TuningResultTable()
) {

    fun addRecord(name: String, parameter: String, error: LabelErrors, optimal: Boolean) {
        results.addRecord(name, parameter, error)
        if (optimal) {
            optimalResults.addRecord(name, parameter, error)
        }
    }

    fun saveOptimalResults(optimalParametersPath: Path) {
        optimalResults.save(optimalParametersPath)
        LOG.info("Saved optimal tuning results to: $optimalParametersPath")
        LOG.info("Optimal tuning error: ${optimalResults.minError()}")
    }

    fun saveTuningErrors(resultsPath: Path) {
        results.save(resultsPath)
        LOG.info("Saved all the tuning results to: $resultsPath")
    }

    companion object {
        private val LOG = LoggerFactory.getLogger(TuningResults::class.java)
    }
}
