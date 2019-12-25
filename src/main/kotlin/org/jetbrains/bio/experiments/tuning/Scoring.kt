package org.jetbrains.bio.experiments.tuning

import org.apache.log4j.Logger
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.containers.LocationsMergingList
import java.nio.file.Path
import java.util.*


/**
 * Simple counter of correct out of all
 */
data class ErrorRate(var total: Int = 0, var correct: Int = 0) {

    fun observe(outcome: Boolean) {
        total++; if (outcome) correct++
    }

    fun combine(other: ErrorRate) {
        total += other.total; correct += other.correct
    }

    fun rate(): Double = 1.0 * correct / total
}


/**
 * Merging map based error [LocationLabel] -> [ErrorRate]
 */
class LabelErrors(val map: MutableMap<LocationLabel, ErrorRate> = TreeMap())
    : MutableMap<LocationLabel, ErrorRate> by map {

    fun error(type: Label? = null) = 1.0 - correct(type) * 1.0 / total(type)

    fun correct(type: Label?) = filter(type)
            .stream()
            .mapToInt { (_, errRate) -> errRate.correct }
            .sum()

    fun total(type: Label?) = filter(type)
            .stream()
            .mapToInt { (_, errRate) -> errRate.total }
            .sum()

    private fun filter(type: Label?) = when (type) {
        null -> map.entries
        else -> map.entries.filter { (ann, _) -> ann.type == type }
    }

    fun combine(other: LabelErrors) {
        other.map.entries.forEach {
            if (it.key !in map) {
                map[it.key] = ErrorRate()
            }
            map[it.key]?.combine(it.value)
        }
    }

    fun observe(label: LocationLabel, outcome: Boolean) {
        var ler = map[label]
        if (ler == null) {
            ler = ErrorRate()
            map[label] = ler
        }
        ler.observe(outcome)
    }
}

class TuningResultTable {

    val names: MutableList<String> = arrayListOf()
    private val parameters: MutableList<String> = arrayListOf()
    private val errors: Array<MutableList<Double>> = Array(Label.values().size + 1) {
        arrayListOf<Double>()
    }

    fun minError() = errors.mapNotNull { it.min() }.min()

    fun addRecord(name: String, parameter: String, labelErrors: LabelErrors) {
        names.add(name)
        parameters.add(parameter)
        Label.values().forEachIndexed { i, type -> errors[i].add(labelErrors.error(type)) }
        errors[Label.values().size].add(labelErrors.error())
    }

    fun print(path: Path) {
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


class TuningResults(private val results: TuningResultTable = TuningResultTable(),
                    private val optimalResults: TuningResultTable = TuningResultTable()) {


    fun addRecord(name: String, parameter: String, error: LabelErrors, optimal: Boolean) {
        results.addRecord(name, parameter, error)
        if (optimal) {
            optimalResults.addRecord(name, parameter, error)
        }
    }

    fun saveOptimalResults(optimalParametersPath: Path) {
        optimalResults.print(optimalParametersPath)
        LOG.info("Saved optimal tuning results to: $optimalParametersPath")
        LOG.info("Optimal tuning error: ${optimalResults.minError()}")
    }

    fun saveTuningErrors(resultsPath: Path) {
        results.print(resultsPath)
        LOG.info("Saved all the tuning results to: $resultsPath")
    }

    companion object {
        private val LOG = Logger.getLogger(TuningResults::class.java)
    }
}

/**
 * Function to compute Peaks [peaks] vs Labels [labels]
 * @return map like structure
 */
fun computeErrors(labels: List<LocationLabel>, peaks: LocationsMergingList): LabelErrors {
    val labelErrors = LabelErrors()
    for (label in labels) {
        labelErrors.observe(label, label.check(peaks))
    }
    return labelErrors
}
