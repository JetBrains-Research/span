package org.jetbrains.bio.tools

import org.apache.log4j.Logger
import org.jetbrains.bio.dataset.ChipSeqTarget
import org.jetbrains.bio.dataset.DataConfig
import org.jetbrains.bio.genome.Genome
import org.jetbrains.bio.util.*
import java.io.BufferedReader
import java.io.IOException
import java.io.InputStreamReader
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.Path
import java.util.*


/**
 * Main entry point to the washu repository scripts
 * See https://github.com/JetBrains-Research/washu
 */
class Washu(private val washuPath: Path = PATH) {

    companion object {
        val PATH: Path
            get() {
                val path = "/mnt/stripe/washu".toPath()
                check(path.exists && path.isDirectory) {
                    "FAILED to find washu repository at: $path"
                }
                return path
            }

        internal val LOG = Logger.getLogger(Washu::class.java)
    }

    /**
     * @return absolute path for given [relativePath]
     */
    private fun script(relativePath: String): String {
        val path = "$washuPath/$relativePath".toPath()
        check(path.exists && path.isReadable) {
            "Missing file $path in $washuPath"
        }
        return path.toAbsolutePath().toString()
    }

    fun runMACS2(genome: Genome,
                 output: Path,
                 fdr: Double,
                 broadPeaks: Boolean = false,
                 extParams: List<String> = listOf()) {
        LOG.info("Working directory: $output")

        val suffix = (if (broadPeaks) {
            "broad"
        } else {
            "fdr"
        }) + "$fdr"

        val additionalMacs2Args = arrayListOf("-B")
        if (broadPeaks) {
            additionalMacs2Args.add("--broad")
            additionalMacs2Args.add("--broad-cutoff")
            additionalMacs2Args.add(fdr.toString())
        } else {
            additionalMacs2Args.add("-q")
            additionalMacs2Args.add(fdr.toString())
        }
        additionalMacs2Args.addAll(extParams)

        val args = arrayListOf(
                genome.build,
                "-", // Don't pass CHROM_SIZES argument to prevent big file creation
                suffix,
                "'${additionalMacs2Args.joinToString(separator = " ")}'",
                output.toString())

        val scriptPath = output / "run_MACS2_$suffix.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("parallel/macs2.sh")} ${args.joinToString(separator = " ")}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }


    fun runSicer(genome: Genome,
                 output: Path,
                 fdr: Double,
                 extParams: List<String> = listOf()) {
        LOG.info("Working directory: $output")

        var args = listOf(output.toString(),
                genome.build,
                genome.chromSizesPath.toString(),
                fdr.toString())

        args += extParams

        val scriptPath = output / "run_Sicer.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("parallel/sicer.sh")} ${args.joinToString(separator = " ")}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runRSeg(genome: Genome,
                output: Path) {
        LOG.info("Working directory: $output")

        val args = listOf(output.toString(),
                genome.build,
                genome.chromSizesPath.toString())


        val scriptPath = output / "run_Rseg.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("parallel/rseg.sh")} ${args.joinToString(separator = " ")}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runManorm(genome: Genome,
                  readsPath: Path,
                  peaksPath: Path,
                  output: Path,
                  name: String) {
        val chromSize = genome.chromSizesPath.toString()
        val csvConfig = createDiffBindConfig(readsPath, peaksPath, output, name)
        val scriptPath = output / "run_manorm_diff.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_manorm.sh")} $name $chromSize $csvConfig")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runDiffReps(genome: Genome,
                    readsPath: Path,
                    output: Path,
                    name: String,
                    broadPeaks: Boolean) {
        val chromSize = genome.chromSizesPath.toString()
        val scriptPath = output / "run_diffreps.sh"

        val args = if (broadPeaks) "--mode block --nsd broad" else "-me nb"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_diffreps.sh")} $name $chromSize $readsPath \"$args\"")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runDiffBind(output: Path,
                    name: String,
                    csvConfig: Path) {
        val scriptPath = output / "run_diffbind_diff.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_diffbind.sh")} $name $csvConfig")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }


    fun runIntersection(output: Path, vararg beds: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$washuPath && " +
                "bash ${script("bed/intersect.sh")} ${beds.map { it.toAbsolutePath() }.joinToString(" ")} > $output") {
            directory(output.parent.toFile())
        }
    }

    fun runUnion(output: Path, vararg beds: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$washuPath && " +
                "bash ${script("bed/union.sh")} ${beds.map { it.toAbsolutePath() }.joinToString(" ")} > $output") {
            directory(output.parent.toFile())
        }
    }

    fun runReads2BW(readsPath: Path, chromSizesPath: Path, output: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$washuPath && " +
                "bash ${script("scripts/reads2bw.sh")} $readsPath $chromSizesPath $output") {
            directory(output.parent.toFile())
        }
    }

    fun runReads2TagsBW(readsPath: Path, fragmentSize: Int, chromSizesPath: Path, output: Path) {
        "bash".exec("-c", "export WASHU_ROOT=$washuPath && " +
                "bash ${script("downstream/signals/bam2tagsbw.sh")} $readsPath $fragmentSize $chromSizesPath $output") {
            directory(output.parent.toFile())
        }
    }

    fun runPoolReads(paths: List<Path>, chromSizesPath: Path, outputBam: Path) {
        val scriptPath = outputBam.parent / "pool_reads.sh"
        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("BAMS=()")
            writer.newLine()
            paths.forEach {
                writer.write("BAM=$(bash ${script("scripts/reads2bam.sh")} $it $chromSizesPath); BAMS+=(\"\$BAM\")")
                writer.newLine()
            }
            writer.write("samtools merge $outputBam \${BAMS[@]}")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(outputBam.toFile())
        }
    }

    fun runDiffMacsPooled(name: String, output: Path) {
        val scriptPath = output / "run_macs_diff_pooled.sh"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_macs_diff_pooled.sh")} $name $output")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }

    }

    fun runChipDiff(name: String, genome: Genome, output: Path, macsPooled: Path) {
        val scriptPath = output / "run_chipdiff_pooled.sh"

        val chromSize = genome.chromSizesPath.toString()

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_chipdiff.sh")} $name $output $chromSize $macsPooled")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runMacsBdgdiff(name: String, output: Path, macsPooled: Path) {
        val scriptPath = output / "run_macs_bdgdiff.sh"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("downstream/diff/run_macs_bdgdiff.sh")} $name $output $macsPooled")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(output.toFile())
        }
    }

    fun runChipDiffPostProcess() {
        "python".exec("$washuPath/downstream/diff/chipseq_diff_postprocess.py")
    }

    fun runRip(bamPath: Path, peakPath: Path) {
        val scriptPath = peakPath.parent / "run_rip.sh"

        scriptPath.bufferedWriter().use { writer ->
            writer.write("export WASHU_ROOT=$washuPath")
            writer.newLine()
            writer.write("bash ${script("scripts/rip.sh")} $bamPath $peakPath")
            writer.newLine()
        }
        "bash".exec(scriptPath) {
            directory(peakPath.parent.toFile())
        }
    }

    /**
     * Nucleotide consensus with regions supported by at least a half of all files.
     */
    fun medianNucleotideConsensus(resultFile: Path, tracks: List<Path>) {
        "bash".exec("-c",
                "bash $washuPath/bed/consensus.sh -p 50 ${tracks.joinToString(separator = " ")} > $resultFile")
    }

    /**
     * Nucleotide consensus with regions supported by at least two files.
     */
    fun weakNucleotideConsensus(resultFile: Path, tracks: List<Path>) {
        "bash".exec("-c",
                "bash $washuPath/bed/consensus.sh -c 2 ${tracks.joinToString(separator = " ")} > $resultFile")
    }
}


fun checkBedtools() {
    val processBuilder = ProcessBuilder("bedtools", "--version")
    try {
        val process = processBuilder.start()
        process.waitFor()
    } catch (e: IOException) {
        e.printStackTrace()
        System.exit(0)
    }
}

fun checkR() {
    val processBuilder = ProcessBuilder("Rscript", "--version")
    try {
        val process = processBuilder.start()
        process.waitFor()
    } catch (e: IOException) {
        e.printStackTrace()
        System.exit(0)
    }
}

fun checkPythonVersion() {
    val processBuilder = ProcessBuilder("python", "--version")
    val process = processBuilder.start()
    val lines = BufferedReader(InputStreamReader(process.errorStream)).use(BufferedReader::readLines)
    process.waitFor()

    val text = lines.joinToString(separator = "\n")
    Washu.LOG.info("Python version: $text")
    val matchResult = "\\d+.\\d+.\\d".toRegex().find(text)!!.value

    if (!matchResult.startsWith("3.")) {
        System.err.println("python version found $matchResult")
        System.err.println("python command should point to Python 3")
        System.exit(1)
    }
}

/**
 * Creates symbolic links in the [basePath] given Chip-Seq [modification]
 */
fun DataConfig.fillData(basePath: Path,
                                                   modification: String?,
                                                   addInput: Boolean = true,
                                                   addFailedTracks: Boolean = false): List<Path> {
    val files = ArrayList<Path>()
    for ((key, section) in tracksMap) {
        for ((replicate, contents) in section.filter { addFailedTracks || !it.second.failedTrack }) {
            val isInput = ChipSeqTarget.isInput(key.dataType)
            val addTrack = if (isInput) addInput else (modification == null || modification == key.dataType)
            if (addTrack) {
                val localName = "${key.cellId}_${
                if (isInput)
                    "input"
                else
                    "${replicate}_${key.dataType}"
                }.bam"
                val localPath = basePath / localName
                if (Files.exists(localPath, LinkOption.NOFOLLOW_LINKS)) {
                    localPath.delete()
                }
                Washu.LOG.info("${contents.path} -> $localPath")
                Files.createSymbolicLink(localPath, contents.path)
                files.add(localPath)
            }
        }
    }
    return files
}

/**
 * Creates symbolic links using [fillData] and launches [body] code.
 */
fun DataConfig.runBatch(path: Path, mark: String?,
                                                   addInput: Boolean = true,
                                                   addFailedTracks: Boolean = false,
                                                   body: (Path) -> Unit) {
    path.checkOrRecalculate(isDirectory = true) {
        val basePath = it.path
        val files = fillData(basePath, mark, addInput = addInput, addFailedTracks = addFailedTracks)
        body(basePath)
        files.forEach { Files.delete(it) }
    }
}

/**
 * Depends on the data layout in [fillData]
 */
fun createDiffBindConfig(readsPath: Path, peaksPath: Path, output: Path,
                         name: String, tissue: String = "CD14", factor: String = "Age"): Path {

    return createDiffBindConfig(readsPath, output, name, tissue, factor) { replicate, read_path ->
        val peaks = (peaksPath.glob("*$replicate*.bed") +
                peaksPath.glob("*${read_path.fileName.stem}*.broadPeak") +
                peaksPath.glob("*${read_path.fileName.stem}*.narrowPeak")
                ).filterNot { "pileup" in it.toString() }

        check(peaks.isNotEmpty()) {
            "No peak files found for $replicate $read_path in $peaksPath"
        }
        check(peaks.size == 1) {
            Washu.LOG.warn("More than 1 peak file found for $replicate $read_path in $peaksPath\n$peaks")
        }
        peaks.single()
    }
}


fun createDiffBindConfig(readsPath: Path, output: Path,
                         name: String, tissue: String = "CD14", factor: String = "Age",
                         peaksProvider: (String, Path) -> Path): Path {
    val csvConfig = output / "$name.csv"
    csvConfig.bufferedWriter().use {
        it.write("SampleID,Tissue,Factor,Condition,Replicate,bamReads,ControlID,bamControl,Peaks,PeakCaller")
        it.newLine()
        readsPath.glob("*.bam")
                .filter { "input" !in it.toString() }
                .sorted()
                .forEach { READ ->
                    val split = READ.fileName.toString().split("_")
                    // Use O or Y as a condition for washu chipseq_diff.sh compatibility
                    val CONDITION = split[0].substring(0, 1)
                    val REPLICATE = split[1]
                    val SAMPLEID = REPLICATE
                    val INPUT = readsPath.glob("*$CONDITION*input*.bam").single()
                    val replicate = READ.fileName.stem.replace("_unique", "")
                    val PEAK = peaksProvider(replicate, READ)
                    it.write("$SAMPLEID,$tissue,$factor,$CONDITION,$REPLICATE,$READ,${CONDITION}_input,$INPUT,$PEAK,bed")
                    it.newLine()
                }
    }
    Washu.LOG.info(csvConfig.bufferedReader().readLines().joinToString("\n"))
    return csvConfig
}
