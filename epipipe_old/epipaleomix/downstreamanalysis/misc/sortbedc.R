args <- commandArgs(trailingOnly = TRUE)
###args <- c('in', 'out')

readdf <- function(f){
   read.table(f)
}

df <- readdf(args[1])
n <- do.call(rbind,strsplit(with(df, sprintf('%s_%s_%s',V1,V2,V5)),'_'))
mode(n) = 'numeric'
df = df[order(n[,3],n[,4],n[,5],n[,1],n[,2]),]
gz1 <- gzfile(args[2], "w")
write.table(df, file=gz1 ,row.names=F,col.names=F,quote=F,sep='\t')
close(gz1)
