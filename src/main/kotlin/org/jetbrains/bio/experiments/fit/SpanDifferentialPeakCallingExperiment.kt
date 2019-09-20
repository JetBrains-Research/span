package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.SpanFitInformation.Companion.chromSizes
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.containers.genomeMap
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.Query
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.span.CoverageScoresQuery
import org.jetbrains.bio.span.Peak
import org.jetbrains.bio.span.getChromosomePeaks
import org.jetbrains.bio.span.scoresDataFrame
import org.jetbrains.bio.statistics.MultiLabels
import org.jetbrains.bio.statistics.hmm.MLConstrainedNBHMM
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.state.ZLHID
import org.jetbrains.bio.util.div
import java.nio.file.Path

/**
 * @author Alexey Dievsky
 * @since 10/04/15
 */
class SpanDifferentialPeakCallingExperiment(
        genomeQuery: GenomeQuery,
        paths1: List<SpanDataPaths>,
        paths2: List<SpanDataPaths>,
        fragment: Fragment,
        binSize: Int,
        unique: Boolean
): SpanModelFitExperiment<MLConstrainedNBHMM, Span1CompareFitInformation, ZLHID>(
    Span1CompareFitInformation.effective(
        genomeQuery,
        paths1, paths2,
        MultiLabels.generate(TRACK1_PREFIX, paths1.size).toList(),
        MultiLabels.generate(TRACK2_PREFIX, paths2.size).toList(),
        fragment, unique, binSize
    ),
    MLConstrainedNBHMM.fitter(paths1.size, paths2.size),
    MLConstrainedNBHMM::class.java,
    ZLHID.values(), NullHypothesis.of(ZLHID.same())
) {

    override val defaultModelPath: Path = experimentPath / "${fitInformation.id}.span"

    fun computeDirectedDifferencePeaks(fdr: Double,
            gap: Int): Pair<List<Peak>, List<Peak>> {
        val map = genomeMap(genomeQuery, parallel = true) { chromosome ->
            results.getChromosomePeaks(chromosome, fdr, gap, dataQuery.apply(chromosome))
        }
        val highLow = arrayListOf<Peak>()
        val lowHigh = arrayListOf<Peak>()
        genomeQuery.get().forEach { chromosome ->
            val states = getStatesDataFrame(chromosome)
            map[chromosome].forEach {
                if (states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.D.name) >
                        states.getAsFloat(it.startOffset / fitInformation.binSize, ZLHID.I.name)) {
                    highLow.add(it)
                } else {
                    lowHigh.add(it)
                }
            }
        }
        return highLow to lowHigh
    }


    companion object {
        const val TRACK1_PREFIX = "track1_"
        const val TRACK2_PREFIX = "track2_"

        /**
         * Creates experiment for model-based comparison of binned coverage tracks for given queries.
         *
         * @return experiment [SpanDifferentialPeakCallingExperiment]
         */
        fun getExperiment(
                genomeQuery: GenomeQuery,
                paths1: List<SpanDataPaths>,
                paths2: List<SpanDataPaths>,
                bin: Int,
                fragment: Fragment,
                unique: Boolean
        ): SpanDifferentialPeakCallingExperiment {
            check(paths1.isNotEmpty() && paths2.isNotEmpty()) { "No data" }
            return SpanDifferentialPeakCallingExperiment(genomeQuery, paths1, paths2, fragment, bin, unique)
        }
    }
}

data class Span1CompareFitInformation(
        override val build: String,
        val data1: List<SpanDataPaths>,
        val data2: List<SpanDataPaths>,
        val labels1: List<String>,
        val labels2: List<String>,
        val fragment: Fragment,
        val unique: Boolean,
        override val binSize: Int,
        override val chromosomesSizes: LinkedHashMap<String, Int>
) : SpanFitInformation {

    override val id: String = reduceIds(
        data1.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                listOf("vs") +
                data2.flatMap { listOfNotNull(it.treatment, it.control) }.map { it.stemGz } +
                listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
    )

    override val dataQuery: Query<Chromosome, DataFrame>
        get() {
            val data = data1 + data2
            val labels = labels1 + labels2
            return object : CachingQuery<Chromosome, DataFrame>() {
                val scores = data.map {
                    CoverageScoresQuery(genomeQuery(), it.treatment, it.control, fragment, binSize, unique)
                }

                override fun getUncached(input: Chromosome): DataFrame {
                    return scores.scoresDataFrame(input, labels.toTypedArray())
                }

                override val id: String
                    get() = reduceIds(scores.zip(labels).flatMap { (s, l) -> listOf(s.id, l) })
            }
        }


    override fun scoresDataFrame(): Map<Chromosome, DataFrame> {
        val gq = genomeQuery()
        val queries1 = data1.map {
            CoverageScoresQuery(gq, it.treatment, it.control, fragment, binSize, unique)
        }
        val queries2 = data2.map {
            CoverageScoresQuery(gq, it.treatment, it.control, fragment, binSize, unique)
        }
        if (queries1.any { !it.ready } || queries2.any { !it.ready }) {
            return emptyMap()
        }
        return gq.get().associateBy({it}) {
            DataFrame.columnBind(
                queries1.scoresDataFrame(it, labels1.toTypedArray()),
                queries2.scoresDataFrame(it, labels2.toTypedArray())
            )
        }
    }

    companion object {
        const val VERSION: Int = 3

        fun effective(
                genomeQuery: GenomeQuery,
                paths1: List<SpanDataPaths>,
                paths2: List<SpanDataPaths>,
                labels1: List<String>,
                labels2: List<String>,
                fragment: Fragment,
                unique: Boolean,
                binSize: Int
        ): Span1CompareFitInformation {
            return Span1CompareFitInformation(
                genomeQuery.build, paths1, paths2, labels1, labels2, fragment, unique, binSize,
                chromSizes(
                    SpanModelFitExperiment.effectiveGenomeQuery(
                        genomeQuery, paths1 + paths2, fragment, unique
                    )
                )
            )
        }
    }
}