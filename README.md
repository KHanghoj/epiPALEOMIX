# epiPALEOMIX

A Fast, Accurate, and Automatic pipeline for generating nucleosome and methylation maps from high throughput sequencing data underlying ancient samples.

## Brief description

*   [Overview](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-Overview)
*   [Installation](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-Installation)
*   [How to run epiPALEOMIX](https://bitbucket.org/khanghoj/epiomix/overview#markdown-header-How-to-run-epiPALEOMIX)
*   [Tutorial](tutorial.html)

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

*   Python 2.7.3-10
*   R v2.15+
*   [SAMTools](http://samtools.sourceforge.net) v0.1.18+ ([Li et al 2009](http://bioinformatics.oxfordjournals.org/content/25/16/2078.long))

*   Python modules:
    *   [pysam](https://github.com/pysam-developers/pysam/) v0.8+

            $ easy_install pip  # Pip is a python package manager.
            $ pip install pysam


## epiPALEOMIX pipeline installation

*   Install all required dependencies as listed above.

* Download and extract the [epiPALEOMIX archive](https://bitbucket.org/khanghoj/epiomix/overview). You can use the command-lines below:

        $ git clone https://khanghoj@bitbucket.org/khanghoj/epiomix.git
        or
        $ wget https://khanghoj@bitbucket.org/khanghoj/epiomix.git

* You can optionally create a symbolic link to add epiPALEOMIX in your executable paths to avoid writing the full path. For example, if ~/bin is in your executable path ($PATH):

        $ cd ~/bin
        $ ln -s -T path/to/epiomix/epaleomix.py epiPALEOMIX

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

For an example run, see [Tutorial](tutorial.html)

## Explaining the output of epiPALEOMIX

Two folder are generated when running the epiPALEOMIX pipeline. Assuming no `--destination` option added two folders are created in the current working directory:

1. "OUT__makefileName_" contains One folder for each BAM-file analyzed containing final output files generated by epiPALEOMIX in _"BAMName\_AnalysisName\_BedName.txt.gz"_ format.
2. A "TEMP__makefile_" folder containing all temporary files generated by epiPALEOMIX. No Final results are located in this folder.

The epiPALEOMIX pipeline produces flat tabulated output files in gzip format. Output files of methylation analyses contain genomic position, counts of deaminated reads, coverage. Nucleosome maps are provided as output files showing the genomic coordinates for each predicted nucleosome, peak read-depth, and nucleosome calling score. Phasogram output files contain the distribution of distances between successive read starts within a window of size defined by users. (4) WriteDepth contains genomic position, read depth, running nucleosome score. The package also implements a read-depth correction procedure for variation in %GC content. When the GC-correction option is enabled, the model created and used in the analyses is provided as an additional output together with the GC-corrected read depth for each genomic position.




### for tutorial
*   No extra R-Packages are required to run epiPALEOMIX itself. However, attached Rscripts in [example](example) requires the following:
    *   [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
    *   [reshape2](https://cran.r-project.org/web/packages/reshape2/index.html)

