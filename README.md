# epiPALEOMIX Pipeline

A Fast, Accurate, and Automatic pipeline for generating nucleosome and methylation maps from high throughput sequencing data underlying ancient samples.


Recent developments in ancient DNA research have shown that not only ancient genomes but ancient epigenomes can be recovered from high-throughput DNA sequence data. Although direct methods, such as bisulfite DNA treatment [(Llamas et al. 2012)][llamas] and target enrichment for methylated epialleles [(Seguin-Orlando et al. 2015)][seguin], have been used, only indirect computational methods exploiting patterns of DNA degradation along the genome have succeeded in reconstructing genome-wide epigenetic maps [(Pedersen et al. 2014][pedersen], [Gokhman et al. 2014)][gokhman]. epiPALEOMIX implements the latter methods in a user-friendly and open-source pipeline, allowing non-experienced users to reconstruct both methylation and nucleosome maps from ancient DNA data. It requires BAM alignment files against reference genomes as input, which can easily be generated using the PALEOMIX pipeline, previously released by our group [(Schubert et al. 2014)][schubert]. epiPALEOMIX has been developed by Kristian Hanghøj in the research group of Dr. Ludovic Orlando, at the Centre for GeoGenetics, Natural History Museum of Denmark, University of Copenhagen, Denmark.


For detailed instructions and an in-depth description of epiPALEOMIX, please refer to the companion Wiki site.

For questions, bug reports, and/or suggestions contact Kristian Hanghøj at k.hanghoej@snm.ku.dk

## Citation

Hanghøj K., Seguin-Orlando A., Schubert M., Madsen H.T., Pedersen J.S., Willerslev E., Orlando L. A Fast, Accurate, and Automatic pipeline for generating nucleosome and methylation maps from high throughput sequencing data underlying ancient samples (2016) Molecular Biology and Evolution, in press.

## References

+ Llamas, B. et al. High-resolution analysis of cytosine methylation in ancient DNA. PLoS One 7, e30226 (2012).
+ Seguin-Orlando, A. et al. Pros and cons of methylation-based enrichment methods for ancient DNA. Sci. Rep. 5, 11826 (2015).
+ Pedersen, J. S. et al. Genome-wide nucleosome map and cytosine methylation levels of an ancient human genome. Genome Res. 24, 454–466 (2014).
+ Gokhman, D. et al. Reconstructing the DNA methylation maps of the Neandertal and the Denisovan. Science 344, 523–527 (2014).
+ Schubert, M. et al. Characterization of ancient and modern genomes by SNP detection and phylogenomic and metagenomic analysis using PALEOMIX. Nat. Protoc. 9, 1056–1082 (2014).


***
[Go to wiki](https://bitbucket.org/khanghoj/epipaleomix/wiki/Home)
***

[llamas]: http://dx.doi.org/10.1371/journal.pone.0030226
[seguin]: http://dx.doi.org/10.1038/srep11826
[pedersen]: http://dx.doi.org/10.1101/gr.163592.113
[gokhman]: http://dx.doi.org/10.1126/science.1250368
[schubert]: http://dx.doi.org/10.1038/nprot.2014.063