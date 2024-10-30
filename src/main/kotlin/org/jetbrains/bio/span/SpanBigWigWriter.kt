package org.jetbrains.bio.span

import org.jetbrains.bio.CompressionType
import org.jetbrains.bio.big.BigWigFile
import org.jetbrains.bio.big.FixedStepSection
import org.jetbrains.bio.big.WigSection
import org.jetbrains.bio.genome.GenomeQuery
import org.jetbrains.bio.genome.Strand
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.util.Progress
import org.jetbrains.bio.util.await
import org.jetbrains.bio.util.time
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import org.slf4j.event.Level
import java.nio.file.Path
import java.util.concurrent.Callable
import java.util.concurrent.atomic.AtomicLong

object SpanBigWigWriter {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun write(
        spanResults: SpanFitResults,
        genomeQuery: GenomeQuery,
        bigWigPath: Path,
        blackListPath: Path? = null
    ) {
        val fitInfo = spanResults.fitInfo
        check(fitInfo is SpanAnalyzeFitInformation) {
            "Cannot create bigwig coverage is only available for analyze command"
        }
        fitInfo.prepareData()
        check(fitInfo.normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        check(fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }) {
            "Coverage information is not available"
        }

        val blackList = if (blackListPath != null) {
            LOG.info("Loading blacklist regions: $blackListPath")
            LocationsMergingList.load(genomeQuery, blackListPath)
        } else null

        val chromosomes = genomeQuery.get()
        val scoresSumProgress = Progress {
            title = "Estimating total scores"
        }.bounded(chromosomes.size.toLong())
        val scoresSum = AtomicLong()
        chromosomes.map { chromosome ->
            Callable {
                val data = fitInfo.dataQuery.apply(chromosome)
                val scores = data.sliceAsInt(data.labels.first())
                val chrSum = scores.sumOf { it.toLong() }
                scoresSum.addAndGet(chrSum)
                scoresSumProgress.report()
            }
        }.await(true)

        val scale = 1e6 / scoresSum.get()
        LOG.debug("Total scores: ${scoresSum.get()}, scale: $scale")


        val collectDataProgress = Progress {
            title = "Creating bigwig file"
        }.bounded(chromosomes.size.toLong())

        val fakeSection = FixedStepSection("chrFake", 100, span = 100)
        val wigSections = Array<WigSection>(chromosomes.size) {
            fakeSection
        }

        val bin = fitInfo.binSize
        try {
            chromosomes.mapIndexed { i, chromosome ->
                Callable {
                    val data = fitInfo.dataQuery.apply(chromosome)
                    val scores = data.sliceAsInt(data.labels.first())
                    if (blackList != null) {
                        blackList[chromosome, Strand.PLUS].forEach { (startOffset, endOffset) ->
                            for (b in startOffset / bin until endOffset / bin) {
                                scores[b] = 0
                            }
                        }
                        blackList[chromosome, Strand.MINUS].forEach { (startOffset, endOffset) ->
                            for (b in startOffset / bin until endOffset / bin) {
                                scores[b] = 0
                            }
                        }
                    }
                    val section = FixedStepSection(
                        chromosome.name, 0,
                        step = bin, span = bin
                    ).apply {
                        scores.forEach { value -> this.add((value * scale).toFloat()) }
                    }
                    wigSections[i] = section
                    collectDataProgress.report()
                }
            }.await(true)
        } finally {
            collectDataProgress.done()
        }

        LOG.time(Level.INFO, "Saving bigwig scores: $bigWigPath")
        {
            BigWigFile.write(
                wigSections.asIterable(),
                chromosomes.map { it.name to it.length },
                bigWigPath,
                // Old compression for compatibility with IGV 2.3.92 browser
                compression = CompressionType.DEFLATE
            ) {}
        }

    }
}