package org.jetbrains.bio.span

import org.apache.commons.csv.CSVFormat
import org.jetbrains.bio.genome.*
import org.jetbrains.bio.genome.containers.LocationsMergingList
import org.jetbrains.bio.genome.coverage.Roughness
import org.jetbrains.bio.genome.format.BedFormat
import org.jetbrains.bio.genome.format.unpack
import org.jetbrains.bio.span.fit.SpanAnalyzeFitInformation
import org.jetbrains.bio.span.fit.SpanConstants.SPAN_DEFAULT_GAP
import org.jetbrains.bio.span.fit.SpanFitResults
import org.jetbrains.bio.span.peaks.ModelToPeaks.estimateCandidatesNumberLenDist
import org.jetbrains.bio.span.peaks.ModelToPeaks.getChromosomeCandidates
import org.jetbrains.bio.span.peaks.Peak
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_BACKGROUND_SENSITIVITY_VARIANTS
import org.jetbrains.bio.span.semisupervised.SpanSemiSupervised.SPAN_GAPS_VARIANTS
import org.jetbrains.bio.util.*
import org.slf4j.Logger
import org.slf4j.LoggerFactory
import java.io.BufferedWriter
import java.io.FileWriter
import java.nio.file.Path

object SpanResultsAnalysis {

    val LOG: Logger = LoggerFactory.getLogger(javaClass)

    fun doDeepAnalysis(
        actualModelPath: Path,
        spanFitResults: SpanFitResults,
        fitInfo: SpanAnalyzeFitInformation,
        genomeQuery: GenomeQuery,
        fdr: Double,
        blacklistPath: Path?,
        peaksList: List<Peak>,
        peaksPath: Path?
    ) {
        check(fitInfo.normalizedCoverageQueries != null) {
            "Please use prepareData before!"
        }
        check(fitInfo.normalizedCoverageQueries!!.all { it.areCachesPresent() }) {
            "Coverage information is not available"
        }
        val infoFile = if (peaksPath != null) "$peaksPath.txt" else null
        infoFile?.toPath()?.deleteIfExists()
        val infoWriter = if (infoFile != null) BufferedWriter(FileWriter(infoFile)) else null

        // Save basic stats to infoFile
        LOG.info("Processing basic info")
        if (infoWriter != null) {
            val aboutModel = spanFitResults.modelInformation(actualModelPath)
            val aboutPeaks = PeaksInfo.compute(
                genomeQuery,
                peaksList.map { it.location }.stream(),
                peaksPath!!.toUri(),
                fitInfo.paths.map { it.treatment }
            )
            infoWriter.write((aboutModel + aboutPeaks).joinToString("\n") {
                    (k, v) -> "${k.name}: ${k.render(v)}"
            } + "\n")
        }

        LOG.info("Analysing tracks roughness...")
        val ncq = fitInfo.normalizedCoverageQueries!!.first()
        val (treatmentScale, controlScale, beta) = ncq.coveragesNormalizedInfo
        val treatmentCoverage = ncq.treatmentReads.coverage()
        val controlCoverage = ncq.controlReads?.coverage()
        val roughness = Roughness.detectAverageRoughness(
            genomeQuery, treatmentCoverage, treatmentScale, controlCoverage, controlScale, beta
        )
        logInfo("Track roughness: ${"%.3f".format(roughness)}", infoWriter)
        val roughnessNoControl = Roughness.detectAverageRoughness(
            genomeQuery, treatmentCoverage, treatmentScale, null, null, 0.0
        )
        logInfo("Track roughness no control: ${"%.3f".format(roughnessNoControl)}", infoWriter)

        if (blacklistPath != null) {
            LOG.info("Computing blacklisted roughness")
            val roughnessBlackListed = Roughness.detectAverageRoughness(
                genomeQuery, treatmentCoverage, 1.0, null, null, 0.0, blacklistPath
            )
            logInfo(
                "Track roughness blacklisted: ${"%.3f".format(roughnessBlackListed)}",
                infoWriter
            )
        }
        infoWriter?.close()

        LOG.info("Analysing candidates characteristics wrt sensitivity and gap...")
        val sensDetailsFile = if (peaksPath != null) "$peaksPath.sensitivity.tsv" else null
        sensDetailsFile?.toPath()?.deleteIfExists()
        val sensDetailsWriter = if (sensDetailsFile != null)
            BufferedWriter(FileWriter(sensDetailsFile))
        else
            null
        logInfo(
            "Sensitivity\tGap\tCandidatesN\tCandidatesAL\tCandidatesAD", sensDetailsWriter, false
        )
        for (bgs in SPAN_BACKGROUND_SENSITIVITY_VARIANTS) {
            for (g in SPAN_GAPS_VARIANTS) {
                val (candidatesN, candidatesAL, candidatesAD) =
                    estimateCandidatesNumberLenDist(genomeQuery, spanFitResults, fdr, bgs, g, null)
                logInfo(
                    "$bgs\t$g\t$candidatesN\t$candidatesAL\t$candidatesAD", sensDetailsWriter, false
                )
            }
        }
        sensDetailsWriter?.close()

        LOG.info("Analysing peaks segmentation wrt sensitivity")
        val bedtoolsPath = "bedtools".toPath()
        if (peaksPath != null) {
            if (!bedtoolsPath.isAccessible()) {
                LOG.warn("bedtools not installed. Cannot create sensitivity track view.")
            } else {
                val blackList =
                    if (blacklistPath != null) LocationsMergingList.load(genomeQuery, blacklistPath) else null
                withTempDirectory(peaksPath.fileName!!.stem) { dir ->
                    LOG.info("Saving sensitivity peaks to $dir")
                    val peaksPaths = arrayListOf<String>()
                    val peaksSens = doubleArrayOf(10.0, 5.0, 1.0, 0.1, 1e-4, 1e-6)
                    for (s in peaksSens) {
                        val candidates = genomeQuery.get()
                            .flatMap { chr ->
                                if (!spanFitResults.fitInfo.containsChromosomeInfo(chr)) {
                                    return@flatMap emptyList<ChromosomeRange>()
                                }
                                getChromosomeCandidates(
                                    spanFitResults,
                                    chr,
                                    fdr,
                                    s,
                                    SPAN_DEFAULT_GAP
                                ).first.map { it.on(chr) }
                            }
                            .filter { blackList == null || !blackList.intersects(it.location) }
                        val path = dir / "peaks_$s.peak"
                        CSVFormat.TDF.print(path.bufferedWriter()).use { printer ->
                            candidates.forEach { r ->
                                printer.printRecord(r.chromosome.name, r.startOffset.toString(), r.endOffset.toString())
                            }
                        }
                        peaksPaths.add(path.toString())
                    }
                    val out = dir / "multiinter.peak"
                    val sensitivityTrack = "$peaksPath.sensitivity.bed"
                    sensitivityTrack.toPath().deleteIfExists()
                    withTempFile("script", ".sh") { script ->
                        script.bufferedWriter().use { w ->
                            w.write("bedtools multiinter -i ${peaksPaths.joinToString(" ")} | awk -v OFS='\\t' '{print $1, $2, $3, $5}' > $out")
                        }
                        "bash".exec(script.toString())
                    }
                    val scoredLocations = arrayListOf<XLocation<String>>()
                    val format = BedFormat.auto(out)
                    format.parse(out) {
                        it.forEach { entry ->
                            val chromosome = genomeQuery[entry.chrom]
                            if (chromosome != null) {
                                val e = entry.unpack(format)
                                scoredLocations.add(XLocation(Location(e.start, e.end, chromosome), e.name))
                            }
                        }
                    }
                    val lls = Array(peaksPaths.size) { LocationsMergingList.builder(genomeQuery) }
                    scoredLocations.sortedBy { it.l }.forEach { r ->
                        val idx = r.x.split(',').minOf { it.toInt() }
                        lls[idx - 1].add(r.l)
                    }
                    val peaks =
                        lls.flatMapIndexed { i, ll -> ll.build().toList().map { XLocation(it, i) } }.sortedBy { it.l }
                    CSVFormat.TDF.print(sensitivityTrack.toPath().bufferedWriter()).use { printer ->
                        peaks.forEach { r ->
                            val score = (1000.0 * (peaksSens.size - r.x) / peaksSens.size).toInt()
                            val s = peaksSens[r.x]
                            printer.printRecord(
                                r.l.chromosome.name,
                                r.l.startOffset * fitInfo.binSize,
                                r.l.endOffset * fitInfo.binSize,
                                "s$s",
                                score
                            )
                        }
                        LOG.info("Sensitivity file saved to $sensitivityTrack")
                    }
                }
            }
        }
    }

    private fun logInfo(msg: String, infoWriter: BufferedWriter?, useLog: Boolean = true) {
        if (useLog)
            SpanCLA.LOG.info(msg)
        else
            println(msg)
        if (infoWriter != null) {
            infoWriter.write(msg)
            infoWriter.newLine()
        }
    }

    class XLocation<T>(val l: Location, val x: T) : LocationAware {
        override val location: Location
            get() = l
    }

}