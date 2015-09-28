## zcat allgeneshg19.bed.gz | cut -f 1,2,3,4,5 |sed 's/chr//g' | sed '/_/ d' | gzip - > allgeneshg19_cleaned.bed.gz
df <- read.table('allgeneshg19_cleaned.bed.gz', comment.char='!', h=T)
df <- df[!with(df, txEnd-txStart)>1500000,]  # remove super long transcripts greater than 1.5 million 
autosomal <- 1:22
df <- df[(df[,2]%in%autosomal),]
print(table(keep <- !duplicated(sprintf('%s_%s',as.character(df[,2]), as.numeric(df[,4])))))
df <- df[keep,]
print(table(keep <- !duplicated(sprintf('%s_%s',as.character(df[,2]), as.numeric(df[,5])))))
df <- df[keep,]
## keep <- print(table(!duplicated(sprintf('%s_%s_%s',as.character(df[,2]), as.numeric(df[,4]),as.numeric(df[,5])))))
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
df <-  df[,c(2,4,5,3,1)]
colnames(df) <- c('chrom', 'start', 'end', 'strand', 'ID')


plusstrand <- df$strand == '+'
minusstrand <- df$strand == '-'

plusdf <- with(df, data.frame(chrom=chrom[plusstrand], start=start[plusstrand]-500,
                              end=start[plusstrand]+1000, GEBObegin=start[plusstrand]+1000,
                              GEBOend=end[plusstrand], strand=strand[plusstrand], ID=ID[plusstrand]))

minusdf <- with(df, data.frame(chrom=chrom[minusstrand], start=end[minusstrand]-1000,
                               end=end[minusstrand]+500, GEBObegin=start[minusstrand],
                               GEBOend=end[minusstrand]-1000,
                               strand=strand[minusstrand],
                               ID=ID[minusstrand]))

df <- rbind(plusdf, minusdf)
df$PROMREG=sprintf('%s_%s_%s',df[,1], df[,2], df[,3])
df$GEBOREG=sprintf('%s_%s_%s',df[,1], df[,4], df[,5])
print(table(keep <- !duplicated(df$PROMREG)))
df = df[keep,]
df <- df[with(df, (GEBObegin-GEBOend)<0),]

write.table(df[,c(1,2,3,6,7)], file='PROM_autosom_wochr.bed',
            row.names=F,col.names=F,quote=F,sep='\t')
write.table(df[,c(1,4,5,6,7)], file='GEBO_autosom_wochr.bed',
            row.names=F,col.names=F,quote=F,sep='\t')

#df$PROMREG=sprintf('%s_%s_%s',df[,1], df[,2], df[,3])
#df$GEBOREG=sprintf('%s_%s_%s',df[,1], df[,4], df[,5])
write.table(df, file='GEBOPROMTOGETHER.INFOFILE',
            row.names=F,col.names=T,quote=F,sep='\t')


cat('remember to sort files:\n sort -k 1,1n -k 2,2n -k 3,3n PROM_autosom_wochr.bed > tmp.file  && mv tmp.file PROM_autosom_wochr.bed\n')


cat('remember to sort files:\n sort -k 1,1n -k 2,2n -k 3,3n GEBO_autosom_wochr.bed > tmp.file  && mv tmp.file GEBO_autosom_wochr.bed\n')
