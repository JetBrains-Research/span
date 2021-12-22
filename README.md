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

**SPAN Peak Analyzer** is a multipurpose peak caller capable of processing a broad range of ChIP-seq, ATAC-seq, and
single-cell ATAC-seq datasets.<br>
In [semi-supervised mode](https://artyomovlab.wustl.edu/aging/tools) it is capable to robustly handle multiple
replicates and noise by leveraging limited manual annotation information.

**Open Access Paper:** https://doi.org/10.1093/bioinformatics/btab376

**Citation:** Shpynov O, Dievskii A, Chernyatchik R, Tsurinov P, Artyomov MN. Semi-supervised peak calling with SPAN and
JBR Genome Browser. Bioinformatics. 2021 May 21.

The Latest release
--------------
Version [0.13.5244](https://github.com/JetBrains-Research/span/releases/tag/0.13.5244) released on Aug 12th, 2020.

Requirements
------------

Download and install [Java 8](http://www.java.com/en/download/).

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

SPAN Command line options
-------------------------

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

`-b, --bin BIN_SIZE`<br>
Peak analysis is performed on read coverage tiled into consequent bins of configurable size. (default: 200)

`-f, --fdr FDR`<br>
Minimum FDR cutoff to call significant regions. (default: 0.05)

`-g, --gap GAP`<br>
Gap size to merge spatially close peaks. Useful for wide histone modifications. (default: 3)

`--labels LABELS`<br>
Labels BED file. Used in semi-supervised peak calling.

`-p, --peaks PEAKS`<br>
Resulting peaks file in ENCODE broadPeak* (BED 6+3) format. If omitted, only the model fitting step is performed.

`-m, --model MODEL`<br>
This option is used to specify SPAN model path, if not provided, model name is composed of input names and other
arguments.

`-w, --workdir PATH`<br>
Path to the working directory (stores coverage and model caches).

`--peaks-type PEAKS_TYPE`<br>
Peaks computation method.<br>
Use 'islands' to merge consequent blocks of enriched bins with relaxed gaps, or 'simple' to merge fdr enriched HMM bins
with gap into peaks (previous). (default: 'islands')

`--threads THREADS`<br>
Configures the parallelism level. SPAN utilizes both multithreading and specialized processor extensions like SSE2, AVX,
etc.

`--ms, --multistarts`<br>
Number of multi-start runs using different model initializations. Use 0 to disable (default: 5)

`--ms-iterations, --msi`<br>
Number of iterations for each multi-start run (default: 2)

`--ms-iterations, --msi`<br>
Maximum number of iterations for EM algorithm. (default: 20)

`--threshold, --tr`<br>
Convergence threshold for EM algorithm, use `--debug` option to see detailed info (default: 1)

`-d, --debug`<br>
Print all the debug information, used for troubleshooting.

`-q, --quiet`<br>
Turn off output.

Example
-------
Step-by-step example with test dataset is available [here](https://github.com/JetBrains-Research/span/wiki).


Galaxy
------

SPAN is available as a tool in the official [ToolShed](https://toolshed.g2.bx.psu.edu/view/jetbrains/span/66b2c9a128ab)
for
[Galaxy](https://galaxyproject.org/). You can ask your Galaxy administrator to install it.

Build from sources
------------------

Clone [bioinf-commons](https://github.com/JetBrains-Research/bioinf-commons) library under the project root.

  ```
  git clone git@github.com:JetBrains-Research/bioinf-commons.git
  ```

Launch the following command line to build SPAN jar:

  ```
  ./gradlew shadowJar
  ```

The SPAN jar file will be generated in the folder `build/libs`.

FAQ
---

* Q: What is the average running time?<br>
  A: SPAN is capable of processing a single ChIP-Seq track in less than 1 hour on an average laptop (MacBook Pro 2015).
* Q: Which operating systems are supported?<br>
  A: SPAN is developed in modern Kotlin programming language and can be executed on any platform supported by java.
* Q: Where did you get this lovely span picture?<br>
  A: From [ascii.co.uk](https://ascii.co.uk), the original author goes by the name jgs.

Errors Reporting
-----------------

Use [GitHub issues](https://github.com/JetBrains-Research/span/issues) to suggest new features or report bugs.

Authors
-------

[JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)