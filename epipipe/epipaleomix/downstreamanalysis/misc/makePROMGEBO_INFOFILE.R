printfunc <- function(df,start,end){
    sprintf('%s_%s_%s', df[,1], df[,start], df[,end])
}

df <- read.table('allgeneshg19.bed.gz', comment.char='!', h=T)
df <- df[!with(df, txEnd-txStart)>1500000,]  # remove super long transcripts greater than 1.5 million 
autosomal <- paste0('chr',1:22)
df <- df[(df$chrom%in%autosomal),]
library(GenomicRanges, quietly = TRUE)  # remove overlapping transcripts
newdf <- GRanges(as.character(df[,2]), IRanges(as.numeric(df[,4]), as.numeric(df[,5])),
                 strand=df$strand, names=df[,1])
newdf <- reduce(newdf)  # takes the longest range
ord <- order(as.character(newdf@seqnames),newdf@ranges)
newdf <- newdf[ord,]
df <- as.data.frame(newdf)[,c(1,2,3,5)] #newdf[(disjointBins(newdf)==1),])
                                        #keep <- findOverlaps(newdf, newdf, select='first')
                                        #df <- df[unique(keep),]
colnames(df) <- c('chrom', 'start', 'end', 'strand')
plusstrand <- df$strand == '+'
minusstrand <- df$strand == '-'

plusdf <- with(df, data.frame(chrom=chrom[plusstrand], start=start[plusstrand]-500,
                              end=start[plusstrand]+1000, GEBObegin=start[plusstrand]+1000,
                              GEBOend=end[plusstrand], strand=strand[plusstrand]))

minusdf <- with(df, data.frame(chrom=chrom[minusstrand], start=end[minusstrand]-1000,
                               end=end[minusstrand]+500, GEBObegin=start[minusstrand],
                               GEBOend=end[minusstrand]-1000, strand=strand[minusstrand]))
df <- rbind(plusdf, minusdf)
keep <- with(df, (GEBObegin-GEBOend)<0)
df <- df[keep,]
df$PROMreg <- printfunc(df, 2,3)
df$GEBOreg <- printfunc(df, 4,5)
print(head(df))
write.table(df, file=sprintf('%s.INFOFILE', 'PROM_GEBO_autosom_wchr'),
            row.names=F,col.names=F,quote=F,sep='\t')
df$chrom = with(df, sapply(strsplit(as.character(chrom),'chr'),'[[',2))
df$PROMreg <- printfunc(df, 2,3)
df$GEBOreg <- printfunc(df, 4,5)
write.table(df, file=sprintf('%s.INFOFILE', 'PROM_GEBO_autosom_wochr'),
            row.names=F,col.names=F,quote=F,sep='\t')


## printfunc(df, 2,3,4,5,'PROM_GEBO_autosom_w0chr')
