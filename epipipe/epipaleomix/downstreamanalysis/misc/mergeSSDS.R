args <- commandArgs(trailingOnly = TRUE)
###args <- c('hmm2.txt.gz', 'hmm6.txt.gz', 'merged.txt.gz')
readdf <- function(f){
    df <- read.table(f)
    df$nampos <- with(df, sprintf('%s_%s_%s',V1,V2,V5))
    df
}

df <- merge(readdf(args[1]), readdf(args[2]), by='nampos', all=TRUE)
n <- do.call(rbind,strsplit(df$nampos,'_'))
final.df <- cbind(n[,1], n[,2],
                  rowSums(df[,c('V3.x','V3.y')],na.rm=T),
                  rowSums(df[,c('V4.x','V4.y')],na.rm=T),
                  paste(n[,3],n[,4],n[,5],sep='_'))

gz1 <- gzfile(args[3], "w")
write.table(final.df, file=gz1 ,row.names=F,col.names=F,quote=F,sep='\t')
close(gz1)
