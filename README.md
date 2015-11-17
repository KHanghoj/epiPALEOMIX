# epiPALEOMIX

A Fast, Accurate, and Automatic pipeline for generating nucleosome and methylation maps from high throughput sequencing data underlying ancient samples.

## Brief description

*   [Overview](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-overview)
*   [Installation](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-installation)
*   [How to run epiPALEOMIX](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-how-to-run-epipaleomix)
*   [Tutorial](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-tutorial.html)
*   [Makefile Documentation](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-documentaion-makefile.html)

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
*   [Bedtools](https://github.com/arq5x/bedtools2/releases)

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

Below is a brief overview of required steps to run epiPALEOMIX after installation as shown in [Installation](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-installation)

## Get help

    $ epiPALEOMIX -h  # Prints extensive options to standard out
    or
    $ epiPALEOMIX help # Prints simple help

## Generate and preparing a makefile

A makefile, provided by the user in [YAML format](http://www.yaml.org/spec/1.2/spec.html), defines all the analyses to be performed in epiPALEOMIX. The makefile contains paths to input files and a specification of paramaters and analyses to be conducted. Full instructions describing the format and structure of epiPALEOMIX Yaml-format makefiles are given in the [Makefile Documentaion](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-documentaion-makefile.html) textfile provided as part of epiPALEOMIX.
epiPALEOMIX requires minimum _one positional argument_ a [yaml-format makefile](http://www.yaml.org/spec/1.2/spec.html). 
To generate a generic makefile `$ epiPALEOMIX makefile > makefile.yaml`, the fill in input paths to a reference genome, bedfile paths, and BAM files and enabled the analyses of interest.

## Running epiPALEOMIX
Prior to starting the analyses using epiPALEOMIX, it is recommendable to check if all executables and input/output files to the node graph are correct with `$ epiPALEOMIX dryrun makefile.yaml`.
Adding `--list-output-files` to the command above, list all outputfiles to standard out.
Assuming you have saved all the required paths in a config file, epiPALEOMIX can be executed using the following simple command:

    $ epiPALEOMIX run makefile.yaml

Results will be saved in the current working directory, in a subdirectory named "OUT_yourmakefile", unless a `--destination` has been given.

## epiPALEOMIX Output

Two folders are generated when running the epiPALEOMIX pipeline.

* "OUT\__makefileName_" folder contains a folder for each BAM-file analyzed with final output files for each analyses conducted in _"BAMName\_AnalysisName\_BedName.txt.gz"_ format.
* "TEMP\__makefileName_" folder containing all temporary files generated by epiPALEOMIX. No Final results are located in this folder.

The epiPALEOMIX pipeline produces flat tabulated output files in gzip format for each BED file.

*   __MethylMap__ output files contain genomic position, counts of deaminated reads, coverage.
*   __Phasogram__ output files contain the distribution of distances between successive read starts within a user defined window.
*   __NucleoMap__ output files contain the genomic coordinates for each predicted nucleosome, peak read-depth, and nucleosome calling score.
*   __WriteDepth__ output files contain genomic position, read depth, running nucleosome score.


The package also implements a read-depth correction procedure for variation in %GC content. If the __GCcorrect__ option is enabled in the makefile, the model created and used in the analyses is provided as an additional output. Analyses applied with GC-correction has a "GCcorr" added to the analyses name.


# Tutorial

## 1\.Investigation the nucleosome occupancy and methylation profile of the CTCF binding sites.
CTCF is an insulator that plays a key role in .... [Fu et al 2004](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000138).
Recent studies of epigenetic regulators such as methylation marks of cytosine at CpG sites and nucleosome occupancy surrounding the CTCF binding sites suggest and out-of-phase pattern between the two regulators. 

No extra R-Packages are required to run epiPALEOMIX itself. However, two packages are required to plot the output in this tutorial.

*   [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
*   [gridBase](https://cran.r-project.org/web/packages/gridBase/index.html)

Both packages can be installed by typing `install.packages(c("ggplot2", "gridBase"))` in R.


We are now ready to test the pipeline on high-throughput sequencing data underlying an ancient sample in the following tutorial

## Downloading and extraction of Tutorial data

Downloading and unpacking the tutorial data (1.3 Gb):

        $ wget FTPFTPFTPFTPFTP
        $ tar xvzf tutorial_epipal.tar.gz
        $ cd tutorial_epipal/

The uncompressed directory contains:

+ __prefix__ folder: Human reference genome (_hs.build37.1.fa_) and a mappability (_GENOME\_51\_50.40000-20000.mappability_) file of assembly hg19
+ __data__ folder: BAM file for testing: *TUTORIAL_CTCF.bam* + .bai index
+ __bedfile__ folder: a Bed file (_CTCF_hg19_wochr.bed_) with genomic coordinates of OCCUPIED CTCF binding sites [Fu et al](http://dx.doi.org/10.1371/journal.pgen.1000138) on hg19 ± 1 kb 
+ __yaml-format makefile__: epiPALEOMIX takes a yaml-format makefile (_example.yaml_) as input containing path to input files and what analyses to be conducted
+ __plottingtools__ folder: A couple of simple R-scripts to plot the output of epiPALEOMIX
+ __creatingplots.sh__: An executable bash script to run the R-scripts in __plottingtools__.

## 2\. Preparing a Yaml makefile

A `example.yaml` makefile for this tutorial is already generated, however, we will go through each section of this makefile. Generic makefiles for epiPALEOMIX can be written to standard out using epipaleomix `$ epiPALEOMIX makefile > new_mkfile.yaml`. For extensive explanation of the generic makefile see [Makefile Documentaion](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-documentaion-makefile.html).

The first line of _example.yaml_ specifies the file type:

    # -*- mode: Yaml; -*-

Second line is a timestamp automatically created by epiPALEOMIX

    # Timestamp: 2015-11-10T11:31:09.819521

### prefix section

    Prefixes:
        --FastaPath: ./prefix/hs.build37.1.fa
        --MappabilityPath: ./prefix/GENOME_51_50.40000-20000.mappability


### bedfile section

    BedFiles:  # AT LEAST ONE BEDFILE IS REQUIRED
        EnabledFilter: True  # True, bed coordinates not overlapping a mappability region with at at least 0.9 in 'UniquenessFilter: 0.9' (user definable)
        CTCF: ./bedfile/CTCF.bed # BEDNAME


### bamfile section

    BamInputs:
        EXAMPLEBAM:  # SET TO A MEANINGFUL BAMNAME 
            BamInfo:
                BamPath: ./data/TUTORIAL_CTCF.BAM
            GCcorrect:
                Enabled: True
            NucleoMap:   # The nucleosome calling tool
                Enabled: True
                Apply_GC_Correction: True  # GCcorrect Enabled must be True
            MethylMap:   # The Methylation mapping tool
                Enabled: True
                --Primes: five  # choose from {both, five, three}. Normally 5' (five) for DS and 5' and 3' (both) for SS
            Phasogram:   # The phasogram analysis 
                Enabled: True
                Apply_GC_Correction: False
                --SubsetPileup: 3
            WriteDepth:
                Enabled: True
                Apply_GC_Correction: True


## 3\. Analyzing CTCF binding sites ±1 kb using epiPALEOMIX.
Next we can run the analyses required the generate the nucleosome positioning, methylation marks, and phasogram.

        $ epiPALEOMiX run example.yaml

Finally, we can plot the data by running __creatingplots.sh__

        $ ./creatingplots.sh

Two PDF files are created, a phasogram plot, (_CTCF\_PHASOGRAM.pdf_) and a plot showing the out-of-phase profile between methylation marks and nucleosome positioning (_CTCF\_Nucleosome\_red\_Methylation_blue.pdf_).



[YAML](http://www.yaml.org/ "YAML official site") is a human-readable format to store data in a hierarchical way. Subcategories are indented with one or more spaces compared to the parent category. Always use **spaces**, and **never use tabs** for indentation. [This section](http://www.yaml.org/spec/1.2/spec.html#Preview) of the documentation will give you a quick overview on YAML.

