[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
[![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Biolabs_Span)/statusIcon.svg)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Biolabs_Span&guest=1)
[![DOI](https://zenodo.org/badge/159559994.svg)](https://doi.org/10.5281/zenodo.15131734)

SPAN Peak Analyzer version 2.0
==============================

```           ,        ,
      __.-'|'-.__.-'|'-.__
    ='=====|========|====='=
    ~_^~-^~~_~^-^~-~~^_~^~^~^
```

**SPAN Peak Analyzer version 2.0** is a universal HMM-based peak caller capable of processing a broad range of ChIP-seq, ATAC-seq,
and single-cell ATAC-seq datasets of different quality.<br> 

Features
--------

* Supports both narrow and broad footprint experiments (ChIP-seq, ATAC-seq, DNAse-seq)
* Produces robust results on datasets of different signal-to-noise ratio, including Ultra-Low-Input ChIP-seq
* Produces highly consistent results in multiple-replicates experiment setup
* Tolerates missing control experiment
* Integrated into the JetBrains Research ChIP-seq
  analysis [pipeline](https://github.com/JetBrains-Research/chipseq-smk-pipeline) from raw reads to visualization and
  peak calling
* Integrated with the [JBR](https://github.com/jetBrains-Research/jbr) Genome Browser, uploaded data model allows for
  interactive visualization and fine-tuning
* _Experimentally_ supports multi-replicated mode and differential peak calling mode
* In [semi-supervised mode](https://artyomovlab.wustl.edu/aging/tools) it is capable to robustly handle multiple
  replicates and noise by leveraging limited manual annotation information.

Latest release
------------------
See [releases](https://github.com/JetBrains-Research/span/releases) section for actual information.

SPAN 2.0 enhancements compared to version 1.0
------------------------------------------

SPAN version 2.0 introduces several key improvements over the original semi-supervised SPAN 1.0, most notably eliminating the need for manual markup annotations.
It now operates in a **fully unsupervised mode** with robust default parameters.<br>
Key changes include:

* **Automated Setup**: SPAN 2.0 no longer requires semi-supervised markup to function. It runs directly with improved default settings.
* **Enhanced Preprocessing**: The data preprocessing pipeline has been redesigned, featuring better control regression and smarter initialization of HMM parameters.
* **Constraint-Driven Model Fitting**: The HMM now includes adaptive constraints for noise floor and signal-to-noise ratio, enhancing robustness across datasets with variable quality.
* **New Peak Detection Framework**: Peak identification now leverages post-model analysis and a unified strategy for extracting peaks from HMM output.
* **Improved Replicates-model**: These enhancements significantly boost performance in replicate-based analyses.
* **Expanded Applicability**: SPAN 2.0 is more effective for diverse data types, including ATAC-seq, CUT&RUN, and CUT&Tag, and supports explicit input format declaration, BigWig signal visualization, and summit calling for fine resolution.

The original SPAN 1.0, which required semi-supervised input, is described in:<br>
<i>Shpynov O, Dievskii A, Chernyatchik R, Tsurinov P, Artyomov MN. Semi-supervised peak calling with SPAN and
JBR Genome Browser. Bioinformatics. 2021 May 21. https://doi.org/10.1093/bioinformatics/btab376</i>

Requirements
------------

Download and install [Java 8+](https://openjdk.org/install/).

Peak calling
------------

To analyze a single (possibly replicated) biological condition use `analyze` command. See details with command:

```bash
$ java -jar span.jar analyze --help
```

The `<output.bed>` file will contain predicted and FDR-controlled peaks in the
ENCODE [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) (BED 6+3) format:

```
<chromosome> <peak start offset> <peak end offset> <peak_name> <score> . <coverage or fold/change> <-log p-value> <-log Q-value>
```

Examples:

* Regular peak calling<br>
  `java -Xmx8G -jar span.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -p Results.peak`
* Semi-supervised peak calling<br>
  `java -Xmx8G -jar span.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -l Labels.bed -p Results.peak`
* Model fitting only<br>
  `java -Xmx8G -jar span.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -m Model.span`

Differential peak calling
-------------------------

_Experimental!_
To compare two (possibly replicated) biological conditions use the `compare`. See help for details:

```bash
$ java -jar span.jar compare --help
```

Command line options
-------------------------

| Parameter                                               | Description                                                                                                                                                                                                                                                          | 
|---------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `-t, --treatment TREATMENT`<br/> **required**           | Treatment file. Supported formats: BAM, BED, or BED.gz file. <br/>If multiple files are provided, they are treated as replicates. <br/>Multiple files should be separated by commas: `-t A,B,C`. <br/>Multiple files are processed as replicates on the model level. |
| `-c, --control CONTROL`                                 | Control file. Multiple files should be separated by commas. <br/>A single control file, or a separate file per each treatment file is required. <br/>Follow the instructions for `-t`, `--treatment`.                                                                |
| `-cs, --chrom.sizes CHROMOSOMES_SIZES`<br/>**required** | Chromosome sizes file for the genome build used in TREATMENT and CONTROL files. <br/>Can be downloaded at [UCSC](https://hgdownload.soe.ucsc.edu/downloads.html).                                                                                                    |
| `-b, --bin BIN_SIZE`                                    | Peak analysis is performed on read coverage tiled into consequent bins of configurable size.                                                                                                                                                                         |
| `-f, --fdr FDR`                                         | False Discovery Rate cutoff to call significant regions.                                                                                                                                                                                                             |
| `-p, --peaks PEAKS`                                     | Resulting peaks file in ENCODE broadPeak* (BED 6+3) format. <br> If omitted, only the model fitting step is performed.                                                                                                                                               |
| `-chr, --chromosomes CHROMOSOMES_LIST`                  | Chromosomes to process, multiple chromosomes should be separated by commas.                                                                                                                                                                                          |
| `--format FORMAT`                                       | Reads file format. Supported: BAM, SAM, CRAM, BED. Text format can be in zip or gzip archive.<br>If not provided, guessed from file extensions.                                                                                                                      |
| `--fragment FRAGMENT`                                   | Fragment size. If provided, reads are shifted appropriately. <br>If not provided, the shift is estimated from the data.<br>`--fragment 0` is recommended for ATAC-Seq data processing.                                                                               |
| `-kd, --keep-duplicates`                                | Keep duplicates. By default, SPAN filters out redundant reads aligned at the same genomic position.<br>Recommended for bulk single cell ATAC-Seq data processing.                                                                                                    |
| `--blacklist BLACKLIST_BED`                             | Blacklisted regions of the genome to be excluded from peak calling results.                                                                                                                                                                                          |
| `--labels LABELS`                                       | Labels BED file. Used in semi-supervised peak calling.                                                                                                                                                                                                               |
| `-m, --model MODEL`                                     | This option is used to specify SPAN model path. Required for further semi-supervised peak calling.                                                                                                                                                                   |
| `-w, --workdir PATH`                                    | Path to the working directory. Used to save coverage and model cache.                                                                                                                                                                                                |
| `--bigwig`                                              | Create beta-control corrected counts per million normalized track.                                                                                                                                                                                                   |
| `--hmm-snr SNR`                                         | Fraction of coverage to estimate and guard signal to noise ratio, `0` to disable constraint check.                                                                                                                                                                   |
| `--hmm-low LOW`                                         | Minimal low state mean threshold, guards against too broad peaks, `0` to disable constraint check.                                                                                                                                                                   |
| `--sensitivity SENSITIVITY`                             | Configures log PEP threshold sensitivity for candidates selection.<br>Automatically estimated from the data, or during semi-supervised peak calling.                                                                                                                 |
| `--gap GAP`                                             | Configures minimal gap between peaks.<br>Generally, not required, but used in semi-supervised peak calling.                                                                                                                                                          |
| `--summits`                                             | Calls summits within peaks.<br>Recommended for ATAC-seq and single-cell ATAC-seq analysis.                                                                                                                                                                           |
| `--f-light LIGHT`                                       | Lightest fragmentation threshold to apply compensation gap.<br>Not available when `gap` is explicitly provided.                                                                                                                                                      |                  
| `--f-hard HARD`                                         | Hardest fragmentation threshold to apply compensation gap.<br>Not available when `gap` is explicitly provided.                                                                                                                                                       |                  
| `--f-speed SPEED`                                       | Fragmentation acceleration threshold to compute gap.<br>Not available when `gap` is explicitly provided.                                                                                                                                                             |                  
| `--clip CLIP_TRESHOLD`                                  | Clip max threshold for fine-tune boundaries according to local signal, `0` to disable.                                                                                                                                                                               |
| `--multiple TEST`                                       | Method applied for multiple hypothesis testing.<br/>`BH` for Benjamini-Hochberg, `BF` for Bonferroni.                                                                                                                                                                |
| `-i, --iterations`                                      | Maximum number of iterations for Expectation Maximisation (EM) algorithm.                                                                                                                                                                                            |
| `--tr, --threshold`                                     | Convergence threshold for EM algorithm, use `--debug` option to see detailed info.                                                                                                                                                                                   |
| `--ext`                                                 | Save extended states information to model file.<br>Required for model visualization in JBR Genome Browser.                                                                                                                                                           |
| `--deep-analysis`                                       | Perform additional track analysis - coverage (roughness) and creates multi-sensitivity bed track.                                                                                                                                                                    |
| `--threads THREADS`                                     | Configure the parallelism level.                                                                                                                                                                                                                                     | |
| `-l, --log LOG`                                         | Path to log file, if not provided, it will be created in working directory.                                                                                                                                                                                          |
| `-d, --debug`                                           | Print debug information, useful for troubleshooting.                                                                                                                                                                                                                 |
| `-q, --quiet`                                           | Turn off standard output.                                                                                                                                                                                                                                            |
| `-kc, --keep-cache`                                     | Keep cache files. By default SPAN creates cache files in working directory and cleans up.                                                                                                                                                                            |

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
  A: SPAN is capable of processing a single ChIP-Seq track in less than 10 minutes on an average laptop.
* Q: Which operating systems are supported?<br>
  A: SPAN is developed in modern [Kotlin](https://kotlinlang.org) programming language and can be executed on any
  platform supported by Java.
* Q: Where did you get this lovely span picture?<br>
  A: From [ascii.co.uk](https://ascii.co.uk), the original author goes by the name jgs.

Errors Reporting
-----------------

Use [GitHub issues](https://github.com/JetBrains-Research/span/issues) to suggest new features or report bugs.

Authors
-------

[JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)
