License [![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
Tests Linux [![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:Biolabs_Span)/statusIcon.svg)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=Biolabs_Span&guest=1)

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
[SPAN](https://artyomovlab.wustl.edu/aging/span.html) is a tool for analyzing ChIP-seq data.
Part of integrated [peak calling](https://artyomovlab.wustl.edu/aging/tools.html) solution.

### Requirements

1. Download and install [Java 8][java8].
2. Download the `<build>.chrom.sizes` chromosome sizes of the organism you want to analyze from the UCSC [website][UCSC].

### Analysis

To analyze a single (possibly replicated) biological condition use `analyze` command. See details with command:

```bash
$ java -jar span.jar analyze --help
```

The `<output.bed>` file will contain predicted and FDR-controlled peaks in the [MACS2](https://github.com/taoliu/MACS) format:
```
<chromosome> <peak start offset> <peak end offset> <peak_name> <score> . <coverage or fold/change> <-log p-value> <-log Q-value>
```

`analyze` reports p- and [Q-values] [q] for the null-hypothesis that a
given bin is not enriched with ChIP-Seq modification. Peaks are formed from a list of truly (in the FDR sense)
**enriched** bins for the analyzed biological condition by thresholding
the Q-value with a threshold `alpha` and merging close peaks using `gap` option to broad ones.
This is equivalent to controlling FDR at level `alpha`.
If control is given it will be used with `fragment_size` to compute coverage for analysis.

### Comparison

The compare two (possibly replicated) biological conditions use the `compare`. See help for details:

```bash
$ java -jar span.jar compare --help
```

The structure of output produced by `compare` is similar to that of `analyze`.
The null-hypotheses however are different. By default `compare` assumes that
there is no difference in **enrichment** for the two biological conditions.

### FAQ

#### What file formats are supported?

`span` supports BED reads and BAM as input.

#### Where did you get this lovely span?

From [ascii.co.uk](http://ascii.co.uk/art/bridges), it seems the original author goes by the name `jgs`.


Authors
-------

* [JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)

[java8]: http://www.java.com/en/download/
[q]: http://en.wikipedia.org/wiki/False_discovery_rate#q-value
[UCSC]: http://hgdownload.cse.ucsc.edu/downloads.html
[releases]: https://github.com/JetBrains-Research/span/releases
[tc]: https://teamcity.jetbrains.com/viewLog.html?buildId=lastSuccessful&buildTypeId=Epigenome_span&tab=artifacts&guest=1
[span_scheme]: https://github.com/JetBrains-Research/span/blob/master/span_scheme.pdf
