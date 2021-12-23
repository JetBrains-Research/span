package org.jetbrains.bio.span.semisupervised

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
class LabelErrors(val map: MutableMap<LocationLabel, ErrorRate> = TreeMap()) :
    MutableMap<LocationLabel, ErrorRate> by map {

    fun error(type: Label? = null) = 1.0 - correct(type) * 1.0 / total(type)

    fun correct(type: Label?) = filter(type).sumOf { (_, errRate) -> errRate.correct }

    fun total(type: Label?) = filter(type).sumOf { (_, errRate) -> errRate.total }

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



