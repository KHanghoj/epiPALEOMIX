#!/usr/bin/env bash

cd /home/krishang/data/methylation/RRBS

filelist="
wgEncodeHaibMethylRrbsOsteoblDukeSitesRep1.bed.gz
wgEncodeHaibMethylRrbsOsteoblDukeSitesRep2.bed.gz
wgEncodeHaibMethylRrbsAg10803UwSitesRep1.bed.gz
wgEncodeHaibMethylRrbsAg10803UwSitesRep2.bed.gz
wgEncodeHaibMethylRrbsBjUwSitesRep1.bed.gz
wgEncodeHaibMethylRrbsBjUwSitesRep2.bed.gz
wgEncodeHaibMethylRrbsMelanoSitesRep1.bed.gz
wgEncodeHaibMethylRrbsBcskin0111002BiochainSitesRep1.bed.gz
wgEncodeHaibMethylRrbsNhdfneoUwSitesRep1.bed.gz
wgEncodeHaibMethylRrbsHsmmDukeSitesRep1.bed.gz
wgEncodeHaibMethylRrbsHsmmtDukeSitesRep1.bed.gz
wgEncodeHaibMethylRrbsGm12878HaibSitesRep1.bed.gz
wgEncodeHaibMethylRrbsGm12878HaibSitesRep1.bed.gz
wgEncodeHaibMethylRrbsGm12891HaibSitesRep1.bed.gz
wgEncodeHaibMethylRrbsGm12891HaibSitesRep2.bed.gz
wgEncodeHaibMethylRrbsSkmcUwSitesRep1.bed.gz
wgEncodeHaibMethylRrbsBcbrainh11058nBiochainSitesRep1.bed.gz
wgEncodeHaibMethylRrbsNhaDukeSitesRep1.bed.gz"

parallel "rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethylRrbs/{} ." ::: $filelist

gunzip *.bed.gz
rename 's/wgEncodeHaibMethylRrbs//' *.bed*
rename 's/SitesRep(\d)/_$1/' *.bed*
rename 's/Duke|Haib//' *.bed*

### hg19 Methyl450

cd /home/krishang/data/methylation/methyl450

filelist="wgEncodeHaibMethyl450BjSitesRep1.bed.gz
wgEncodeHaibMethyl450NhaSitesRep1.bed.gz
wgEncodeHaibMethyl450SkmcSitesRep1.bed.gz
wgEncodeHaibMethyl450Gm12878SitesRep1.bed.gz
wgEncodeHaibMethyl450Gm12891SitesRep1.bed.gz
wgEncodeHaibMethyl450Gm12892SitesRep1.bed.gz
wgEncodeHaibMethyl450Ag10803SitesRep1.bed.gz
"

parallel "rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeHaibMethyl450/{} ." ::: $filelist

rename 's/wgEncodeHaibMethyl450//' *.bed.gz
rename 's/SitesRep(\d)/_$1/' *.bed.gz
#rename 's/Duke|Haib//' *.bed.gz
gunzip *.bed.gz


