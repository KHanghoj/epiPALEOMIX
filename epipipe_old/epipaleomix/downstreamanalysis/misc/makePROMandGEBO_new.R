require(data.table)
df <- fread('zcat allgeneshg19withnm1.bed.gz',data.table=F)
print(table(keep <- !(df$hg19.kgXref.refseq%in%"")))
df <- df[keep,]
df <- df[!with(df, hg19.knownGene.txEnd-hg19.knownGene.txStart)>1500000,]
                                        # remove super long transcripts greater than 1.5 million 
print(table(keep <- !grepl('_',df$hg19.knownGene.om)))
df <- df[keep,]
autosomal <- 1:22
df <- df[(df[,2]%in%autosomal),]
print(table(keep <- !duplicated(sprintf('%s_%s',as.character(df[,2]), as.numeric(df[,4])))))
df <- df[keep,]
print(table(keep <- !duplicated(sprintf('%s_%s',as.character(df[,2]), as.numeric(df[,5])))))
df <- df[keep,]


## print(table(keep <- !duplicated(sprintf('%s_%s_%s',as.character(df[,2]), as.numeric(df[,4]),as.numeric(df[,5])))))
## df <- df[keep,]

## library(GenomicRanges, quietly = TRUE)  # remove overlapping transcripts
## newdf <- GRanges(as.character(df[,2]), IRanges(as.numeric(df[,4]), as.numeric(df[,5])),
##                  strand=df$strand, names=df[,1])
## newdf <- reduce(newdf)  # takes the longest range
## ord <- order(as.character(newdf@seqnames),newdf@ranges)
## newdf <- newdf[ord,]
## df <- as.data.frame(newdf)[,c(1,2,3,5,4)] #newdf[(disjointBins(newdf)==1),])
##                                         #keep <- findOverlaps(newdf, newdf, select='first')
##                                         #df <- df[unique(keep),]
cat('\ncheck the 8th column is correct in inputfile\n')
df <-  df[,c(2,4,5,3,8)]
colnames(df) <- c('chrom', 'start', 'end', 'strand', 'ID')


plusstrand <- df$strand == '+'
minusstrand <- df$strand == '-'

plusdf <- with(df, data.frame(chrom=chrom[plusstrand], start=start[plusstrand]-400,
                              end=start[plusstrand]+1000, GEBObegin=start[plusstrand]+2000,
                              GEBOend=end[plusstrand], strand=strand[plusstrand], ID=ID[plusstrand]))

minusdf <- with(df, data.frame(chrom=chrom[minusstrand], start=end[minusstrand]-1000,
                               end=end[minusstrand]+400, GEBObegin=start[minusstrand],
                               GEBOend=end[minusstrand]-2000,
                               strand=strand[minusstrand],
                               ID=ID[minusstrand]))

## plusdf <- with(df, data.frame(chrom=chrom[plusstrand], start=start[plusstrand]-400,
##                               end=start[plusstrand]+1000, GEBObegin=start[plusstrand]+1000,
##                               GEBOend=end[plusstrand], strand=strand[plusstrand], ID=ID[plusstrand]))

## minusdf <- with(df, data.frame(chrom=chrom[minusstrand], start=end[minusstrand]-1000,
##                                end=end[minusstrand]+400, GEBObegin=start[minusstrand],
##                                GEBOend=end[minusstrand]-1000,
##                                strand=strand[minusstrand],
##                                ID=ID[minusstrand]))

df <- rbind(plusdf, minusdf)
df$PROMREG=sprintf('%s_%s_%s',df[,1], df[,2], df[,3])
df$GEBOREG=sprintf('%s_%s_%s',df[,1], df[,4], df[,5])
print(table(keep <- !duplicated(df$PROMREG)))
df <- df[keep,]
print(table(keep <- with(df, (GEBObegin-GEBOend)<0)))
df <- df[keep,]

write.table(df[,c(1,2,3,6,7)], file='PROM_autosom_wochr.bed',
             row.names=F,col.names=F,quote=F,sep='\t')
write.table(df[,c(1,4,5,6,7)], file='GEBO_autosom_wochr.bed',
             row.names=F,col.names=F,quote=F,sep='\t')

write.table(df, file='GEBOPROMTOGETHER.INFOFILE',
            row.names=F,col.names=T,quote=F,sep='\t')


cat('remember to sort files:\n sort -k 1,1n -k 2,2n -k 3,3n PROM_autosom_wochr.bed > tmp.file  && mv tmp.file PROM_autosom_wochr.bed\n')


cat('remember to sort files:\n sort -k 1,1n -k 2,2n -k 3,3n GEBO_autosom_wochr.bed > tmp.file  && mv tmp.file GEBO_autosom_wochr.bed\n')

exit('no')
## this is a special short GEBO. only 1400 bp from tss+2000 to tss+3400
plusstrand <- df$strand == '+'
minusstrand <- df$strand == '-'


plusdf <- with(df, data.frame(chrom=chrom[plusstrand], start=start[plusstrand]-400,
                              end=start[plusstrand]+1000, GEBObegin=start[plusstrand]+2000,
                              GEBOend=start[plusstrand]+3400,
                              strand=strand[plusstrand], ID=ID[plusstrand]))

minusdf <- with(df, data.frame(chrom=chrom[minusstrand], start=end[minusstrand]-1000,
                               end=end[minusstrand]+400,
                               GEBObegin=end[minusstrand]-3400,
                               GEBOend=end[minusstrand]-2000,
                               strand=strand[minusstrand],
                               ID=ID[minusstrand]))

df <- rbind(plusdf, minusdf)
df$PROMREG=sprintf('%s_%s_%s',df[,1], df[,2], df[,3])
df$GEBOREG=sprintf('%s_%s_%s',df[,1], df[,4], df[,5])
print(table(keep <- !duplicated(df$PROMREG)))
df <- df[keep,]

normalprom = read.table('PROM_autosom_wochr.bed')
normalprom$PROMREG=sprintf('%s_%s_%s',normalprom[,1], normalprom[,2], normalprom[,3])
keep <- intersect(normalprom$PROMREG,df$PROMREG)
df <-  df[(df$PROMREG%in%keep),]

write.table(df[,c(1,4,5,6,7)], file='GEBOSHORT_autosom_wochr.bed',
            row.names=F,col.names=F,quote=F,sep='\t')


write.table(df, file='GEBOPROMTOGETHER_GEBOSHORT.INFOFILE',
            row.names=F,col.names=T,quote=F,sep='\t')


cat('remember to sort files:\n sort -k 1,1n -k 2,2n -k 3,3n GEBOSHORT_autosom_wochr.bed > tmp.file  && mv tmp.file GEBOSHORT_autosom_wochr.bed\n')
