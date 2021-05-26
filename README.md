[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
[![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Biolabs_Span)/statusIcon.svg)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Biolabs_Span&guest=1)

SPAN Peak Analyzer
==================

```
+----------------------------------+
|SPAN Semi-supervised Peak Analyzer|
+----------------------------|/----+
           ,        ,
      __.-'|'-.__.-'|'-.__
    ='=====|========|====='=
    ~_^~-^~~_~^-^~-~~^_~^~^~^
```

**SPAN Peak Analyzer** is a semi-supervised multipurpose peak caller capable of processing a broad range of ChIP-seq,
ATAC-seq, and single-cell ATAC-seq datasets that robustly handles multiple replicates and noise by leveraging limited
manual annotation information.\
Part of semi-supervised [peak calling](https://artyomovlab.wustl.edu/aging/tools) solution.

**Open Access Paper:** https://doi.org/10.1093/bioinformatics/btab376

**Citation:** Shpynov O, Dievskii A, Chernyatchik R, Tsurinov P, Artyomov MN. Semi-supervised peak calling with SPAN and JBR Genome Browser. Bioinformatics. 2021 May 21.

The Latest release
--------------
Version [0.13.5244](https://github.com/JetBrains-Research/span/releases/tag/0.13.5244) released on Aug 12th, 2020.

Requirements
------------

1. Download and install [Java 8][java8].
2. Download the `<build>.chrom.sizes` chromosome sizes of the organism you want to analyze from the UCSC [website][UCSC]
   .

Peak calling
------------

To analyze a single (possibly replicated) biological condition use `analyze` command. See details with command:

```bash
$ java -jar span.jar analyze --help
```

The `<output.bed>` file will contain predicted and FDR-controlled peaks in the
ENCODE [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) (BED 6+3) format
(like [MACS2](https://github.com/taoliu/MACS)):

```
<chromosome> <peak start offset> <peak end offset> <peak_name> <score> . <coverage or fold/change> <-log p-value> <-log Q-value>
```

`analyze` reports p- and [Q-values] [q] for the null-hypothesis that a given bin is not enriched with ChIP-Seq
modification. Peaks are formed from a list of truly (in the FDR sense)
**enriched** bins for the analyzed biological condition by thresholding the Q-value with a threshold `alpha` and merging
close peaks using `gap` option to broad ones. This is equivalent to controlling FDR at level `alpha`. If control is
given it will be used with `fragment_size` to compute coverage for analysis.

Examples:

* Regular peak calling<br>
  `java -Xmx4G -jar span.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -p Results.peak`
* Semi-supervised peak calling<br>
  `java -Xmx4G -jar span.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -l Labels.bed -p Results.peak`
* Model fitting only<br>
  `java -Xmx4G -jar span.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes`

Differential peak calling
-------------------------

The compare two (possibly replicated) biological conditions use the `compare`. See help for details:

```bash
$ java -jar span.jar compare --help
```

The structure of output produced by `compare` is similar to that of `analyze`. The null-hypotheses however are
different. By default `compare` assumes that there is no difference in **enrichment** for the two biological conditions.

SPAN Command line options
-------------------------

`-b, --bin BIN_SIZE`<br>
Peak analysis is performed on read coverage tiled into consequent bins of configurable size. Default value is 200bp,
approximately the length of one nucleosome.

`-t, --treatment TREATMENT`<br>
**Required**. ChIP-seq treatment file. Supported formats: BAM, BED, BED.gz or bigWig file. If multiple files are
provided, they are treated as replicates. Multiple files should be separated by commas: `-t A,B,C`. Multiple files are
processed as replicates on the model level.

`-c, --control CONTROL`<br>
Control file. Multiple files should be separated by commas. A single control file, or a separate file per each treatment
file is required. Follow the instructions for `-t`, `--treatment` TREATMENT.

`-cs, --chrom.sizes CHROMOSOMES_SIZES`<br>
**Required**. Chromosome sizes file for the genome build used in TREATMENT and CONTROL files. Can be downloaded
at [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html).

`--fragment FRAGMENT`<br>
Fragment size. If provided, reads are shifted appropriately. If not provided, the shift is estimated from the data.
`--fragment 0` argument is necessary for ATAC-Seq data processing.

`-k, --keep-dup`<br>
Keep duplicates. By default, SPAN filters out redundant reads aligned at the same genomic position.
`--keep-dup` argument is necessary for single cell ATAC-Seq data processing.

`-m, --model MODEL`<br>
This option is used to specify SPAN model path, if not provided, model name is composed of input names and other
arguments.

`-p, --peaks PEAKS`<br>
Resulting peaks file in ENCODE broadPeak* (BED 6+3) format. If omitted, only the model fitting step is performed.

`-f, --fdr FDR`<br>
Minimum FDR cutoff to call significant regions, default value is 1.0E-6. SPAN reports p- and q- values for the null
hypothesis that a given bin is not enriched with a histone modification. Peaks are formed from a list of truly (in the
FDR sense) enriched bins for the analyzed biological condition by thresholding the Q-value with a cutoff `FDR`
and merging spatially close peaks using `GAP` option to broad ones. This is equivalent to controlling FDR. q-values are
calculated from p-values using the Benjamini-Hochberg procedure.

`-g, --gap GAP`<br>
Gap size to merge spatially close peaks. Useful for wide histone modifications. <br>
Default value is 5, i.e. peaks separated by 5 *`BIN` distance or less are merged.

`--labels LABELS`<br>
Labels BED file. Used in semi-supervised peak calling.

`-d, --debug`<br>
Print all the debug information, used for troubleshooting.

`-q, --quiet`<br>
Turn off output.

`-w, --workdir PATH`<br>
Path to the working directory (stores coverage and model caches).

`--threads THREADS`<br>
Configures the parallelism level. SPAN utilizes both multithreading and specialized processor extensions like SSE2, AVX,
etc. Parallel computations were performed using an open-source library [viktor]() for parallel matrix computations in
Kotlin programming language.

`--ms, --multistarts`<br>
Number of multistart runs using different model initializations. Use 0 to disable (default: 5)

`--ms-iterations, --msi`<br>
Number of iterations for each multistart run (default: 5)

`--threshold, --tr`<br>
Convergence threshold for EM algorithm, use `--debug` option to see detailed info (default: 0.1)

Example
-------
Step-by-step example with test dataset is available [here](https://github.com/JetBrains-Research/span/wiki).


Study Cases
-----------
As a benchmark we applied the SPAN peak calling approach to public conventional ChIP-seq datasets as well as to a ULI
ChIP-seq dataset. CD14+ classical monocytes tracks available from the ENCODE database were a natural choice for a
conventional ChIP-seq dataset.<br>
SPAN produced high quality peak calling in all of these cases,
see [report](https://artyomovlab.wustl.edu/aging/study_cases.html).

Galaxy
------

SPAN is available as a tool in the official [ToolShed](https://toolshed.g2.bx.psu.edu/view/jetbrains/span/66b2c9a128ab)
for
[Galaxy](https://galaxyproject.org/). You can ask your Galaxy administrator to install it.

FAQ
---

* Q: What is the average running time?<br>
  A: SPAN is capable of processing a single ChIP-Seq track in less than 1 hour on an average laptop (MacBook Pro 2015).
* Q: Which operating systems are supported?<br>
  A: SPAN is developed in modern Kotlin programming language and can be executed on any platform supported by java.
* Q: Where did you get this lovely span picture?<br>
  A: From [ascii.co.uk](https://ascii.co.uk), the original author goes by
  the name jgs.

Errors Reporting
-----------------

Use this [Issues Tracker](https://github.com/JetBrains-Research/span/issues) to suggest new features or report bugs.

Authors
-------

[JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)

[java8]: http://www.java.com/en/download/

[q]: http://en.wikipedia.org/wiki/False_discovery_rate#q-value

[UCSC]: http://hgdownload.cse.ucsc.edu/downloads.html

[releases]: https://github.com/JetBrains-Research/span/releases

[tc]: https://teamcity.jetbrains.com/viewLog.html?buildId=lastSuccessful&buildTypeId=Epigenome_span&tab=artifacts&guest=1

[span_scheme]: https://github.com/JetBrains-Research/span/blob/master/span_scheme.pdf

