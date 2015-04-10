args <- commandArgs(trailingOnly = TRUE)
mappa_path = args[1]
UNIQUENESS = as.double(args[2])
bed_path = args[3]
dest_path = args[4]
print(mappa_path)
print(UNIQUENESS)
# vec = c()
# for (arg in args[3:length(args)]){
# 	vec = c(vec, arg)
# }
library(GenomicRanges,quietly = TRUE, warn.conflicts=FALSE, verbose=FALSE)
library(data.table,quietly = TRUE, warn.conflicts=FALSE, verbose=FALSE)
# UNIQUENESS = 0.9
# bed_path = 'testhg19ctcf.bed'
# mappa_path = 'chrom22_51_50.40000-20000.mappability'
# mappa_path = 'GENOME_51_50.40000-20000.mappability'

mappa <- data.frame(fread(mappa_path))
mappa$V1 = paste0('chr', mappa$V1)
lastcol <- length(mappa[1,])
table(mappakeep <- mappa[,lastcol]>UNIQUENESS)
mappa <- mappa[mappakeep,]
GR_mappa = GRanges(seqnames = mappa[,1], IRanges(start=mappa[,2], end=mappa[,3]))
bed <- data.frame(fread(bed_path))
GR_bed = GRanges(seqnames = bed[,1], IRanges(start=bed[,2], end=bed[,3]))
hits = countOverlaps(GR_bed, GR_mappa)
keep <-hits>0
bed_new = bed[keep,]
col = dim(bed_new)[2]
write.table(bed_new, file=dest_path,row.names=F,col.names=F,quote=F,sep='\t')
# write(bed_new, stdout(), ncolumns=col)


# test.fun = function(bed_path) {
# 	bed <- data.frame(fread(bed_path))
# 	GR_bed = GRanges(seqnames = bed[,1], IRanges(start=bed[,2], end=bed[,3]))
# 	hits = countOverlaps(GR_bed, GR_mappa)
# 	# write(table(keep <-hits>0),stderr())
# 	keep <- hits>0
# 	bed_new = bed[keep,]
# 	write(bed_new, stdout(), ncolumn=4)
# 	# write.table(bed_new)
# }
# sapply(vec, test.fun)