package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.experiments.fit.FitInformation.Companion.chromSizes
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.sequence.CpGContent
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.query.stemGz
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.mixture.PoissonRegressionMixture
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.util.div
import org.jetbrains.bio.viktor.asF64Array
import java.nio.file.Path

class Span2PeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths: SpanDataPaths,
        mapabilityPath: Path?,
        fragment: Fragment,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>,
        unique: Boolean = true,
        fixedModelPath: Path? = null
) : SpanModelFitExperiment<Model, Span2FitInformation, State>(
    genomeQuery to createDataQuery(genomeQuery, paths, mapabilityPath, fragment, binSize, unique),
    Span2FitInformation(
        genomeQuery.build, listOf(paths),
        listOfNotNull("y", paths.control?.let { "input" }),
        mapabilityPath, fragment, unique, binSize, chromSizes(genomeQuery)
    ),
    modelFitter, modelClass,
    states, nullHypothesis,
    fixedModelPath
) {

    override val modelPath: Path get() = fixedModelPath ?: experimentPath / "$id.span2"
    override val id = reduceIds(
        listOfNotNull(paths.treatment, paths.control, mapabilityPath).map { it.stemGz } +
                listOfNotNull(fragment.nullableInt, binSize).map { it.toString() }
    )

    companion object {

        fun binnedCoverageAsInt(chr: Chromosome, coverage: Coverage, binSize: Int): IntArray {
            val len = (chr.length - 1) / binSize + 1
            val res = IntArray(len)
            for (i in 0 until len - 1) {
                res[i] = coverage.getBothStrandsCoverage(
                    ChromosomeRange(i * binSize, (i + 1) * binSize, chr)
                )
            }
            res[len - 1] = coverage.getBothStrandsCoverage(
                ChromosomeRange((len - 1) * binSize, chr.length, chr)
            )
            return res
        }

        fun binnedCoverageAsDouble(chr: Chromosome, coverage: Coverage, binSize: Int): DoubleArray {
            val len = (chr.length - 1) / binSize + 1
            val res = DoubleArray(len)
            for (i in 0 until len - 1) {
                res[i] = coverage.getBothStrandsCoverage(
                    ChromosomeRange(i * binSize, (i + 1) * binSize, chr)
                ).toDouble()
            }
            res[len - 1] = coverage.getBothStrandsCoverage(
                ChromosomeRange((len - 1) * binSize, chr.length, chr)
            ).toDouble()
            return res
        }

        fun binnedMapability(chr: Chromosome, mapabilityPath: Path, binSize: Int): DoubleArray {
            val bwFile = BigWigFile.read(mapabilityPath)

            val len = (chr.length - 1) / binSize + 1

            if (!bwFile.chromosomes.containsValue(chr.name)) {
                // the chromosome isn't present in the bigWig file, use mean genome mapability for all bins
                val meanMapability = bwFile.totalSummary.sum / bwFile.totalSummary.count
                return DoubleArray(len) { meanMapability }
            }

            val mapSummary = bwFile.summarize(chr.name, 0, (len - 1) * binSize, numBins = len - 1)
            val res = DoubleArray(len)
            for (i in 0 until len - 1) {
                res[i] = mapSummary[i].sum / binSize
            }
            res[len - 1] = bwFile.summarize(chr.name, (len - 1) * binSize, chr.length).single().sum /
                    (chr.length - (len - 1) * binSize)
            return res
        }

        fun createDataQuery(
                genomeQuery: GenomeQuery,
                paths: SpanDataPaths,
                mapabilityPath: Path?,
                fragment: Fragment,
                binSize: Int,
                unique: Boolean
        ) = object : CachingQuery<Chromosome, DataFrame>() {

            private val treatmentCoverage =
                    ReadsQuery(genomeQuery, paths.treatment, unique, fragment, logFragmentSize = false)
            private val controlCoverage = paths.control?.let {
                ReadsQuery(genomeQuery, it, unique, fragment, logFragmentSize = false)
            }

            override val id: String
                get() = reduceIds(listOfNotNull(treatmentCoverage.id, controlCoverage?.id))

            override fun getUncached(input: Chromosome): DataFrame {
                val y = binnedCoverageAsInt(input, treatmentCoverage.get(), binSize)
                val control = controlCoverage?.let { binnedCoverageAsDouble(input, it.get(), binSize) }
                val gc = CpGContent.binnedMeanCG(input, binSize)
                val gc2 = (gc.asF64Array() * gc.asF64Array()).data
                val mapability = mapabilityPath?.let { binnedMapability(input, it, binSize) }
                var df = DataFrame().with("y", y)
                df = df.with("GC", gc).with("GC2", gc2)
                if (control != null) df = df.with("input", control)
                if (mapability != null) df = df.with("mapability", mapability)
                return df
            }
        }

        fun getExperiment(
                genomeQuery: GenomeQuery,
                data: List<SpanDataPaths>,
                mapabilityPath: Path?,
                fragment: Fragment,
                binSize: Int,
                unique: Boolean,
                fixedModelPath: Path?
        ): Span2PeakCallingExperiment<PoissonRegressionMixture, ZLH> {
            check(data.size == 1) { "Poisson regression mixture currently accepts a single data track." }
            return Span2PeakCallingExperiment(
                genomeQuery, data.single(), mapabilityPath,
                fragment, binSize,
                PoissonRegressionMixture.fitter(), PoissonRegressionMixture::class.java,
                ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                unique,
                fixedModelPath
            )
        }
    }
}

data class Span2FitInformation(
        override val build: String,
        override val data: List<SpanDataPaths>,
        override val labels: List<String>,
        val mapabilityPath: Path?,
        override val fragment: Fragment,
        override val unique: Boolean,
        override val binSize: Int,
        override val chromosomesSizes: LinkedHashMap<String, Int>
): SpanAnalyzeFitInformation {
    constructor(
            genomeQuery: GenomeQuery,
            data: List<SpanDataPaths>,
            labels: List<String>,
            mapabilityPath: Path?,
            fragment: Fragment,
            unique: Boolean,
            binSize: Int
    ): this(
        genomeQuery.build, data, labels, mapabilityPath, fragment, unique, binSize, chromSizes(genomeQuery)
    )
}