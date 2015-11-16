# epiPALEOMIX

A Fast, Accurate, and Automatic pipeline for generating nucleosome and methylation maps from high throughput sequencing data underlying ancient samples.

## Brief description

*   [Overview](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-overview)
*   [Installation](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-installation)
*   [How to run epiPALEOMIX](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-how-to-run-epipaleomix)
*   [Tutorial](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-tutorial.html)

# Overview
The epiPALEOMIX pipeline is a free and open-source pipeline for generating two epigenetic regulators of expression data from ancient samples using high-throughput seqeuncing data by benefiting from natural _post-mortem_ degradation processes affecting the DNA. More specifically, the pipeline automates generation of genome-wide methylation marks in CpG context, nucleosome calling and phasogram analyses as shown in [_pedersen et al 2008_](http://dx.doi.org/10.1101/gr.163592.113)

The pipeline takes BAM-format ([Binary form of SAM format](http://bioinformatics.oxfordjournals.org/content/25/16/2078)) files as input containing the aligned reads underlying an ancient sample.

## What it does
epiPALEOMIX: an open-source pipeline for ancient epigenomic analyses
This package is written in [Python 2.7](https://www.python.org/) and builds on the implementation of the [PALEOMIX pipeline](https://github.com/MikkelSchubert/paleomix/wiki/Overview), with makefiles in yaml format (http://yaml.org/). The epiPALEOMIX package handles any type of molecular tools used to prepare aDNA data (including USER-treatment of aDNA extracts and Phusion amplification of aDNA libraries). BAM alignment files([Binary form of SAM format files](http://bioinformatics.oxfordjournals.org/content/25/16/2078)), bed coordinates for genomic regions of interest and a reference genome in fasta format are taken as input files. Mappability files can be optionally supplied to restrict nucleosome calling and pileup output to uniquely mappable regions of the genome.

## Citation

Hanghøj K., ... , Orlando L. _Title_ 

# Installation

## epiPALEOMIX requirements

*   [Python 2.7.3-10](https://www.python.org/downloads/)
*   [R v2.15+](https://www.r-project.org/)
*   [SAMTools](http://samtools.sourceforge.net) v0.1.18+ ([Li et al 2009](http://bioinformatics.oxfordjournals.org/content/25/16/2078.long))

*   Python modules:
    *   [pysam](https://github.com/pysam-developers/pysam/) v0.8+

            $ easy_install pip  # Pip is a python package manager.
            $ pip install pysam

## epiPALEOMIX pipeline installation

* Install all required dependencies listed above.

* Download and extract the [epiPALEOMIX archive](https://bitbucket.org/khanghoj/epiomix/overview). You can use the command-lines below:

        $ mkdir ~/install
        $ cd ~/install
        $ git clone https://khanghoj@bitbucket.org/khanghoj/epiomix.git

* You can optionally create a symbolic link to add epiPALEOMIX in your executable paths to avoid writing the full path. For example, if `~/bin` is in your executable path ($PATH):

        $ cd ~/bin
        $ ln -s -T ~/install/epiomix/epaleomix.py epiPALEOMIX

## Setting default optional arguments

* It is recommended to `--write-config-file` prior to running epiPALEOMIX to get set some default parameters such as temporary folder, number of cores used by default, and warning levels. The epiPALEOMIX.ini file will be written to `~/.pypeline/epiPALEOMIX.ini` and can be edited with any text file editor prior to epiPALEOMIX execution.

        $ epiPALEOMIX --max-threads 20  --write-config-file

* In case a default argument needs to be modified, simply overrule it in the command line, as shown below, or change the `~/.pypeline/epiPALEOMIX.ini` file

        $ epiPALEOMIX run makefile.yaml --max-threads 4

* If a config file is created, the configuration file will be automatically parsed every time epiPALEOMIX is executed, thus the command would be:

        $ epiPALEOMIX run makefile.yaml


# How to run epiPALEOMIX

Below is a brief overview of required steps to run epiPALEOMIX after installation as shown in [Installation](installation.html)

## Get help

    $ epiPALEOMIX -h
    or
    $ epiPALEOMIX help

## Preparing a makefile

epiPALEOMIX requires minimum _one positional argument_ a [yaml-format makefile](http://www.yaml.org/spec/1.2/spec.html). For detailed makefile documentation for epiPALEOMIX see [dingdong](dingdong.html)

## Running epiPALEOMIX

Assuming you have saved all the required paths in a config file, epiPALEOMIX can be executed using the following simple command:

    $ epiPALEOMIX run makefile.yaml

Results will be saved in the current working directory, in a subdirectory named "OUT_yourmakefile", unless you used the "`--destination`" option.

## epiPALEOMIX Output

Two folder are generated when running the epiPALEOMIX pipeline. Assuming no `--destination` option added two folders are created in the current working directory:

1. "OUT\__makefileName_" contains One folder for each BAM-file analyzed containing final output files generated by epiPALEOMIX in _"BAMName\_AnalysisName\_BedName.txt.gz"_ format.
2. A "TEMP\__makefileNAme_" folder containing all temporary files generated by epiPALEOMIX. No Final results are located in this folder.

The epiPALEOMIX pipeline produces flat tabulated output files in gzip format for each BED file.

    * Methylation analysis contain genomic position, counts of deaminated reads, coverage.
    * Phasogram output files contain the distribution of distances between successive read starts within a user defined window.
    * Nucleosome maps are provided as output files showing the genomic coordinates for each predicted nucleosome, peak read-depth, and nucleosome calling score.
    * WriteDepth contains genomic position, read depth, running nucleosome score.
    * The package also implements a read-depth correction procedure for variation in %GC content. If the GC-correction option is enabled, the model created and used in the analyses is provided as an additional output together with the GC-corrected read depth for each genomic position for every sampled.


For an example run, see [Tutorial](tutorial.html)

### for tutorial
*   No extra R-Packages are required to run epiPALEOMIX itself. However, attached Rscripts in [example](example) requires the following:
    *   [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
    *   [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)


# Tutorial

## 1\.Investigation the nucleosome occupancy and methylation profile of the CTCF binding sites.
CTCF is an insulator that plays a key role in .... [Fu et al 2004](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000138).
Recent studies of epigenetic regulators such as methylation marks of cytosine at CpG sites and nucleosome occupancy surrounding the CTCF binding sites suggest and out-of-phase pattern between the two regulators. 

We are now ready to test the pipeline on high-throughput sequencing data underlying an ancient sample in the following tutorial

## Testing epiPALEOMX pipeline

## Downloading and extraction of Tutorial data

Downloading and unpacking the tutorial data (1.3 Gb).


Download and uncompress the data:

        $ wget FTPFTPFTPFTPFTP
        $ tar xvzf tutorial_epipal.tar.gz
        $ cd tutorial_epipal/

The uncompressed directory contains:

1. Prefix folder:
   * Human reference genome, a mappability file of hg19
       * _hs.build37.1.fa_ (Human Reference Genome (HG19/GRCh37)) + .fai index
       * _GENOME_51_50.40000-20000.mappability_, a mappability file of the genome.

2. Data folder:
   * A BAM file for testing
       * _TUTORIAL_CTCF.bam_

3. Bedfile:
   * 
        
1. Reference genome in fasta format (hg19) + .fai index
   * _hs.build37.1.fa_
2. A BAM file for testing + .bai index
   * _TUTORIAL_CTCF.bam_
3. a Bed file with genomic coordinates of CTCF binding sites regions on hg19
   * _CTCF_hg19_nooverlap_wochr.bed_
4. yaml-format makefile
   * _example.yaml_
5. A mappability file of the genome. see [Extras](extras.html) for generation mappability scores.
   * _GENOME_51_50.40000-20000.mappability_
6. R scripts for plotting the output of epiPALEOMIX
   * _XYZ_

## 2\. Preparing a Yaml makefile

A single simple makefile, provided by the user in YAML format, defines all the analyses to be performed in epiPALEOMIX. The makefile contains paths to input files and a specification of paramaters and analyses to be conducted. Full instructions describing the format and structure of epiPALEOMIX Yaml-format makefiles are given in the documentation.makefile textfile provided as part of epiPALEOMIX.

A `example.yaml` makefile for this tutorial is already generated, however, we will go through each section of the makefile. Generic makefiles for epiPALEOMIX can be written to standard out using epipaleomix `$ epiPALEOMIX makefile > new_mkfile.yaml`.

The first line of _example.yaml_ specifies the file type:

    # -*- mode: Yaml; -*-

### Prefix section

### Bedfile section

### Bamfile section

#### Analyses to run


Next we can run the analyses required the generate the nucleosome positioning, methylation marks, and phasogram.


        $ epiPALEOMiX run example.yaml

Finally, we can plot the data using the .R script. Note that ggplot2 is required to run the .R script.

        $ Rscript ????.R

A plot.pdf is generated with the out-of-phase profile of methylation marks and nucleosome prediction occupancy.



[YAML](http://www.yaml.org/ "YAML official site") is a human-readable format to store data in a hierarchical way. Subcategories are indented with one or more spaces compared to the parent category. Always use **spaces**, and **never use tabs** for indentation. [This section](http://www.yaml.org/spec/1.2/spec.html#Preview) of the documentation will give you a quick overview on YAML.

