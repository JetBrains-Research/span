package org.jetbrains.bio.experiments.fit

import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.coverage.Coverage
import org.jetbrains.bio.coverage.Fragment
import org.jetbrains.bio.dataframe.DataFrame
import org.jetbrains.bio.genome.Chromosome
import org.jetbrains.bio.genome.ChromosomeRange
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.sequence.TwoBitSequence
import org.jetbrains.bio.query.CachingQuery
import org.jetbrains.bio.query.ReadsQuery
import org.jetbrains.bio.query.reduceIds
import org.jetbrains.bio.statistics.ClassificationModel
import org.jetbrains.bio.statistics.Fitter
import org.jetbrains.bio.statistics.Preprocessed
import org.jetbrains.bio.statistics.hypothesis.NullHypothesis
import org.jetbrains.bio.statistics.mixture.ZeroPoissonMixture
import org.jetbrains.bio.statistics.state.ZLH
import org.jetbrains.bio.viktor.asF64Array
import java.nio.file.Path

class Span2PeakCallingExperiment<Model : ClassificationModel, State : Any>(
        genomeQuery: GenomeQuery,
        paths: SpanPathsToData,
        fragment: Fragment,
        binSize: Int,
        modelFitter: Fitter<Model>,
        modelClass: Class<Model>,
        states: Array<State>,
        nullHypothesis: NullHypothesis<State>,
        unique: Boolean = true,
        fixedModelPath: Path? = null
) : SpanModelFitExperiment<Model, State>(
    genomeQuery to createDataQuery(genomeQuery, paths, fragment, binSize, unique),
    listOf(paths), listOfNotNull("y", paths.pathInput?.let { "input" }),
    fragment, binSize,
    modelFitter, modelClass,
    states, nullHypothesis,
    unique,
    fixedModelPath
) {

    companion object {

        fun binnedCoverageAsInt(chr: Chromosome, coverage: Coverage, binSize: Int): IntArray {
            val len = (chr.length - 1) / binSize + 1
            val cover = IntArray(len)
            for (i in 0 until len - 1) {
                cover[i] = coverage.getBothStrandsCoverage(ChromosomeRange(i * binSize, (i + 1) * binSize, chr))
            }
            cover[len - 1] = coverage.getBothStrandsCoverage(ChromosomeRange((len-1) * binSize, chr.length, chr))
            return cover
        }

        fun binnedCoverageAsDouble(chr: Chromosome, coverage: Coverage, binSize: Int): DoubleArray {
            val len = (chr.length - 1) / binSize + 1
            val cover = DoubleArray(len)
            for (i in 0 until len - 1) {
                cover[i] = coverage
                        .getBothStrandsCoverage(ChromosomeRange(i * binSize, (i + 1) * binSize, chr))
                        .toDouble()
            }
            cover[len - 1] = coverage
                    .getBothStrandsCoverage(ChromosomeRange((len-1) * binSize, chr.length, chr))
                    .toDouble()
            return cover
        }

        fun meanGC(chr: Chromosome, binSize: Int): DoubleArray {
            val len = (chr.length - 1) / binSize + 1
            val seq: TwoBitSequence = chr.sequence
            val GCcontent = DoubleArray(len)
            for (i in 0 until len - 1) {
                GCcontent[i] = seq.substring(i*binSize, (i + 1)*binSize).count { it == 'c' || it == 'g' }.toDouble()/binSize
            }
            GCcontent[len - 1] = seq
                    .substring((len-1)*binSize, seq.length)
                    .count { it == 'c'|| it == 'g' }
                    .toDouble()/( seq.length - (len-1)*binSize)
            return GCcontent
        }

        fun binnedMapability(chr: Chromosome, mapabilityPath: Path, binSize: Int): DoubleArray {
            if (BigWigFile.read(mapabilityPath).chromosomes.containsValue(chr.name)) {
                val mapSummary = BigWigFile
                        .read(mapabilityPath)
                        .summarize(chr.name, 0, chr.length - chr.length % binSize, numBins = (chr.length - 1) / binSize)
                val result = DoubleArray(mapSummary.size + 1) {
                    if (it < mapSummary.size) mapSummary[it].sum / binSize else 1.0
                }
                result[mapSummary.size] = BigWigFile
                        .read(mapabilityPath)
                        .summarize(chr.name, chr.length - chr.length % binSize, 0)[0].sum / chr.length % binSize
                return result
            }
            val meanMappability = BigWigFile.read(mapabilityPath).totalSummary.sum/ BigWigFile.read(mapabilityPath).totalSummary.count
            return DoubleArray(chr.length) {meanMappability}
        }

        fun createDataQuery(
                genomeQuery: GenomeQuery,
                paths: SpanPathsToData,
                fragment: Fragment,
                binSize: Int,
                unique: Boolean
        ) = object : CachingQuery<Chromosome, DataFrame>() {

            private val treatmentCoverage =
                    ReadsQuery(genomeQuery, paths.pathTreatment, unique, fragment, logFragmentSize = false)
            private val controlCoverage = paths.pathInput?.let {
                ReadsQuery(genomeQuery, it, unique, fragment, logFragmentSize = false)
            }

            override val id: String
                get() = reduceIds(listOfNotNull(treatmentCoverage.id, controlCoverage?.id))

            override fun getUncached(input: Chromosome): DataFrame {
                val y = binnedCoverageAsInt(input, treatmentCoverage.get(), binSize)
                val control = controlCoverage?.let { binnedCoverageAsDouble(input, it.get(), binSize) }
                val gc = meanGC(input, binSize)
                val gc2 = (gc.asF64Array() * gc.asF64Array()).data
                val mapability = paths.pathMappability?.let { binnedMapability(input, it, binSize) }
                var df = DataFrame().with("y", y)
                if (control != null) df = df.with("input", control)
                df = df.with("GC", gc).with("GC2", gc2)
                if (mapability != null) df = df.with("mapability", mapability)
                return df
            }
        }

        private fun semanticCheck(fitter: Fitter<ZeroPoissonMixture>): Fitter<ZeroPoissonMixture> {
            return object : Fitter<ZeroPoissonMixture> by fitter {

                override fun fit(
                        preprocessed: Preprocessed<DataFrame>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): ZeroPoissonMixture =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary()
                        }

                override fun fit(
                        preprocessed: List<Preprocessed<DataFrame>>,
                        title: String,
                        threshold: Double,
                        maxIter: Int,
                        attempt: Int
                ): ZeroPoissonMixture =
                        fitter.fit(preprocessed, title, threshold, maxIter, attempt).apply {
                            flipStatesIfNecessary()
                        }
            }
        }

        fun getExperiment(
                genomeQuery: GenomeQuery,
                data: List<SpanPathsToData>,
                fragment: Fragment,
                binSize: Int,
                unique: Boolean,
                fixedModelPath: Path?
        ): Span2PeakCallingExperiment<ZeroPoissonMixture, ZLH> {
            check(data.size == 1) { "Poisson regression mixture currently accepts a single data track." }
            return Span2PeakCallingExperiment(
                genomeQuery, data.single(),
                fragment, binSize,
                semanticCheck(ZeroPoissonMixture.fitter()), ZeroPoissonMixture::class.java,
                ZLH.values(), NullHypothesis.of(ZLH.Z, ZLH.L),
                unique,
                fixedModelPath
            )
        }
    }
}