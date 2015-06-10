printfunc <- function(df,start,end, prefixname){
    write.table(df[,c(1,start,end,6)], file=sprintf('%s.bed', prefixname),
                row.names=F,col.names=F,quote=F,sep='\t')
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
printfunc(df, 2,3,'PROM_autosom_wchr')
printfunc(df, 4,5,'GEBO_autosom_wchr')
df$chrom = with(df, sapply(strsplit(as.character(chrom),'chr'),'[[',2))
printfunc(df, 2,3,'PROM_autosom_wochr')
printfunc(df, 4,5,'GEBO_autosom_wochr')


## keep <- (minusdf$GEBObegin-minusdf$GEBOend)<0
## minusdf <- minusdf[keep,]
## keep <- (plusdf$GEBObegin-plusdf$GEBOend)<0
## plusdf <- plusdf[keep,]

## printfunc(plusdf, 2,3,'PROM_plusstrand_autosom_wchr')
## printfunc(plusdf, 4,5,'GEBO_plusstrand_autosom_wchr')
## printfunc(minusdf, 2,3,'PROM_minusstrand_autosom_wchr')
## printfunc(minusdf, 4,5,'GEBO_minusstrand_autosom_wchr')


## plusdf[,1] = sapply(strsplit(as.character(plusdf[,1]),'chr'),'[[',2)
## minusdf[,1] = sapply(strsplit(as.character(minusdf[,1]),'chr'),'[[',2)

## printfunc(plusdf, 2,3,'PROM_plusstrand_autosom_wochr')
## printfunc(plusdf, 4,5,'GEBO_plusstrand_autosom_wochr')
## printfunc(minusdf, 2,3,'PROM_minusstrand_autosom_wochr')
## printfunc(minusdf, 4,5,'GEBO_minusstrand_autosom_wochr')
