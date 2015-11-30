# epiPALEOMIX Pipeline

A Fast, Accurate, and Automatic pipeline for generating nucleosome and methylation maps from high throughput sequencing data underlying ancient samples.


Recent developments in ancient DNA research have shown that not only ancient genomes but ancient epigenomes can be recovered from high-throughput DNA sequence data. Although direct methods, such as bisulfite DNA treatment [(Llamas et al. 2012)][llamas] and target enrichment for methylated epialleles [(Seguin-Orlando et al. 2015)][seguin], have been used, only indirect computational methods exploiting patterns of DNA degradation along the genome have succeeded in reconstructing genome-wide epigenetic maps [(Pedersen et al. 2014][pedersen], [Gokhman et al. 2014)][gokhman]. epiPALEOMIX implements the latter methods in a user-friendly and open-source pipeline, allowing non-experienced users to reconstruct both methylation and nucleosome maps from ancient DNA data. It requires BAM alignment files against reference genomes as input, which can easily be generated using the PALEOMIX pipeline, previously released by our group [(Schubert et al. 2013)][schubert]. epiPALEOMIX has been developed by Kristian Hanghøj in the research group of Dr. Ludovic Orlando, at the Centre for GeoGenetics, Natural History Museum of Denmark, University of Copenhagen, Denmark.



For detailed instructions and an in-depth description of epiPALEOMIX, please refer to the companion Wiki site.

For questions, bug reports, and/or suggestions contact Kristian Hanghøj at k.hanghoej@snm.ku.dk

## Citation

Hanghøj K., ... , Orlando L. _Title_ 


## References

+ Llamas, Bastien, Michelle L. Holland, Kefei Chen, Jennifer E. Cropley, Alan Cooper, and Catherine M. Suter. 2012. “High-Resolution Analysis of Cytosine Methylation in Ancient DNA.” PloS One 7 (1): e30226.
+ Seguin-Orlando, Andaine, Cristina Gamba, Clio Der Sarkissian, Luca Ermini, Guillaume Louvel, Eugenia Boulygina, Alexey Sokolov, et al. 2015. “Pros and Cons of Methylation-Based Enrichment Methods for Ancient DNA.” Scientific Reports 5 (July): 11826.
+ Pedersen, Jakob Skou, Eivind Valen, Amhed M. Vargas Velazquez, Brian J. Parker, Morten Rasmussen, Stinus Lindgreen, Berit Lilje, et al. 2014. “Genome-Wide Nucleosome Map and Cytosine Methylation Levels of an Ancient Human Genome.” Genome Research 24 (3): 454–66.
+ Gokhman, David, Eitan Lavi, Kay Prüfer, Mario F. Fraga, José A. Riancho, Janet Kelso, Svante Pääbo, Eran Meshorer, and Liran Carmel. 2014. “Reconstructing the DNA Methylation Maps of the Neandertal and the Denisovan.” Science 344 (6183): 523–27.
+ Schubert, Mikkel, Luca Ermini, Clio Der Sarkissian, Hákon Jónsson, Aurélien Ginolhac, Robert Schaefer, Michael D. Martin, et al. 2014. “Characterization of Ancient and Modern Genomes by SNP Detection and Phylogenomic and Metagenomic Analysis Using PALEOMIX.” Nature Protocols 9 (5): 1056–82.


***
[Go to wiki](https://bitbucket.org/khanghoj/epipaleomix/wiki/Home)
***

[llamas]: http://dx.doi.org/10.1371/journal.pone.0030226
[seguin]: http://dx.doi.org/10.1038/srep11826
[pedersen]: http://dx.doi.org/10.1101/gr.163592.113
[gokhman]: http://dx.doi.org/10.1126/science.1250368
[schubert]: http://dx.doi.org/10.1038/nprot.2014.063

