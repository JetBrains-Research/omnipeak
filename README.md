[![JetBrains Research](https://jb.gg/badges/research.svg)](https://confluence.jetbrains.com/display/ALL/JetBrains+on+GitHub)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg)](https://opensource.org/licenses/MIT)
[![tests](http://teamcity.jetbrains.com/app/rest/builds/buildType:(id:BioLabs_Omnipeak)/statusIcon.svg)](http://teamcity.jetbrains.com/viewType.html?buildTypeId=BioLabs_Omnipeak&guest=1)

OmniPeak
========

**OmniPeak** is a universal HMM-based peak caller capable of processing a broad range of ChIP-seq, ATAC-seq,
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

Requirements
------------

Download and install [Java 21+](https://openjdk.org/install/).

Peak calling
------------

To analyze a single (possibly replicated) biological condition use `analyze` command. See details with command:

```bash
$ java -jar omnipeak.jar analyze --help
```

The `<output.bed>` file will contain predicted and FDR-controlled peaks in the
ENCODE [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13) (BED 6+3) format:

```
<chromosome> <peak start offset> <peak end offset> <peak_name> <score> . <coverage or fold/change> <-log p-value> <-log Q-value>
```

Examples on Java 21:

* Regular peak calling<br>
  `java --add-modules=jdk.incubator.vector -Xmx8G -jar omnipeak.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -p Results.peak`
* Model fitting only<br>
  `java --add-modules=jdk.incubator.vector -Xmx8G -jar omnipeak.jar analyze -t ChIP.bam -c Control.bam --cs Chrom.sizes -m Model.op`

Differential peak calling
-------------------------

_Experimental!_
To compare two (possibly replicated) biological conditions use the `compare`. See help for details:

```bash
$ java --add-modules=jdk.incubator.vector -Xmx8G -jar omnipeak.jar compare --help
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
| `-kd, --keep-duplicates`                                | Keep duplicates. By default, OmniPeak filters out redundant reads aligned at the same genomic position.<br>Recommended for bulk single cell ATAC-Seq data processing.                                                                                                |
| `--blacklist BLACKLIST_BED`                             | Blacklisted regions of the genome to be excluded from peak calling results.                                                                                                                                                                                          |
| `-m, --model MODEL`                                     | This option is used to specify OmniPeak model path.                                                                                                                                                                                                                  |
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
| `-kc, --keep-cache`                                     | Keep cache files. By default OmniPeak creates cache files in working directory and cleans up.                                                                                                                                                                        |

Build from sources
------------------

Clone [bioinf-commons](https://github.com/JetBrains-Research/bioinf-commons) library under the project root.

  ```
  git clone git@github.com:JetBrains-Research/bioinf-commons.git
  ```

Launch the following command line to build OmniPeak jar:

  ```
  ./gradlew shadowJar
  ```

The OmniPeak jar file will be generated in the folder `build/libs`.

Errors Reporting
-----------------

Use [GitHub issues](https://github.com/JetBrains-Research/omnipeak/issues) to suggest new features or report bugs.

Authors
-------

[JetBrains Research BioLabs](https://research.jetbrains.org/groups/biolabs)
