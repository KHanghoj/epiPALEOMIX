library(GenomicRanges)
library(data.table)
UNIQUENESS = 0.9
bed_path = 'testhg19ctcf.bed'
mappa_path = 'chrom22_51_50.40000-20000.mappability'
mappa_path = 'GENOME_51_50.40000-20000.mappability'

mappa <- data.frame(fread(mappa_path))
mappa$V1 = paste0('chr', mappa$V1)
lastcol <- length(mappa[1,])
table(mappakeep <- mappa[,lastcol]>UNIQUENESS)
mappa <- mappa[mappakeep,]

GR_mappa = GRanges(seqnames = mappa[,1], IRanges(start=mappa[,2], end=mappa[,3]))


test.fun = function(bed_path) {
	bed <- data.frame(fread(bed_path))
	GR_bed = GRanges(seqnames = bed[,1], IRanges(start=bed[,2], end=bed[,3]))
	hits = countOverlaps(GR_bed, GR_mappa)
	write(table(keep <-hits>0),stderr())
	bed_new = bed[keep,]
	bed_new
	# write.table(bed_new)
}

abc = test.fun(bed_path)