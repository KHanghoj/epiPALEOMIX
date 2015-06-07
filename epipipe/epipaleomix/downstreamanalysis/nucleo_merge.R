args <- commandArgs(trailingOnly = TRUE)
smooth_val <- as.integer(args[1])
IN_path = args[2]
BEDTYPE = args[3]
OUT_path = args[4]

#  Rscript nucleo_merge.R 50 ~/Desktop/results/allresults/RISESAMPLEOUTPUT/ CTCF ~/Desktop/hmm.pdf
source('~/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/extendsmooth.R')
require(ggplot2) #; require(TSA); require(gridExtra)

read.data.frames <-  function(file, windowsize){
    df <- read.table(file)
    df$bedstart <- sapply(strsplit(as.character(df[,6]), "_"), "[[", 2)
    df$rela_pos <- as.integer(df[,2])-as.integer(df$bedstart)
    df <- df[df[,5]>=1,]
    df <- data.frame(table(df$rela_pos))
    colnames(df) <- c('x','count')
    df$x <- as.numeric(as.character(df$x))
    df$x <- as.integer(max(df$x)/2) - df$x
    df$smooth <- extendsmooth(df$count, windowsize)
    df$plotname <- makename(file)
    as.data.frame(df)
}

plotting <- function(df){
    p <- ggplot(df,aes(x,smooth))+
        geom_line()+ facet_wrap(~plotname,ncol=3,scales='free_y')+
            theme(legend.position="none")+
                labs(title=choosename(BEDTYPE),
                     y='Predicted Nucleosome Counts',
                     x='Relative position')
    p
}

files = list.files(IN_path, pattern='NucleoMap',recursive=T,full.names=T)
files = files[grepl(BEDTYPE, files)]
print(files)
dfs <-  do.call(rbind, parallel::mclapply(files, read.data.frames, windowsize=smooth_val, mc.cores=3))
pdf(OUT_path)
plotting(dfs)
dev.off()
