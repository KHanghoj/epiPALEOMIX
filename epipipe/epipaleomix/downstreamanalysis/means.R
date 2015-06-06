args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]
library(data.table)
require(ggplot2)
source('~/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/extendsmooth.R')
#### Rscript ~/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/means.R ~/Desktop/results/allresults/METHYLSAMPLEOUTPUT/  ~/Desktop/test.pdf

makenames <-  function(f_in, column){
    nam = unlist(strsplit(as.character(f_in),'/'))
    nam = nam[length(nam)]
    as.character(unlist(strsplit(nam,'_'))[column])
}

read.data.frames <-  function(f_in){
    df <- read.table(f_in)
    #df <- data.frame(fread(sprintf('zcat %s', f_in), skip=1, colClasses=list(character=c(1,5))))
    dfs <-  split(df, df[,5])
    df <- data.frame(do.call(rbind, lapply(dfs, f)))
    colnames(df) = c('val', 'totalnumber')
    df$region <- rownames(df)
    df$regionname <- unlist(strsplit(makenames(f_in, 3), '\\.'))[1]
    df$samplename <- makenames(f_in, 1)
#    print(table(keep <- df$totalnumber>=10))
#    df <- df[keep, ] # need a minimum of 10 observations
    print(table(keep <- df$val<quantile(df$val,probs=0.8) & df$totalnumber>=10))
    df[keep, ]  # keeps the values below the 80 % quantile
}

plotting.old <- function(df){
    p <- ggplot(df,aes(x=regionname, y=val, col=regionname, fill=regionname)) +
        geom_bar(position='dodge', stat='identity') +
            facet_wrap(~samplename,ncol=3,scales='free_y') +
                labs(title='', y='MS') + theme(legend.position="none")
    p
}

plotting.newest <- function(df){
    p <- ggplot(df,aes(x=regionname, y=val, col=regionname, fill=regionname)) +
        stat_density(aes(ymax = ..density..,  ymin = -..density.., fill = regionname, color = regionname),
                     geom = "ribbon", position = "identity") +
                         facet_grid(samplename ~ regionname)+coord_flip() + # facet_wrap(~samplename,ncol=3,scales='free_y') +
                             theme(legend.position="none") # + labs(title='', y='MS')
    p
}
plotting <- function(df){
    p <- ggplot(df,aes(x=val)) +
        stat_density(aes(ymax = ..density..,  ymin = -..density.., fill = regionname, color = regionname),
                     geom = "ribbon", position = "identity") +
                         facet_grid(samplename ~ regionname)+coord_flip() + # facet_wrap(~samplename,ncol=3,scales='free_y') +
                             theme(legend.position="none") # + labs(title='', y='MS')
    p
}

plotting.newer <- function(df){
    p <- ggplot(df,aes(x=regionname, y=val, col=regionname, fill=regionname)) +
        geom_boxplot() + facet_wrap(~samplename,ncol=3,scales='free_y') +
            labs(title='', y='MS') + theme(legend.position="none")
    p
}

f <- function(df){
    c(sum(df$V3)/sum(df$V4), sum(df$V4))
}

files <- list.files(IN_path, pattern='MethylMap',recursive=T,full.names=T)
files <- files[grepl('third|PROM|GCSKEWED',files)]
files <- files[grepl('Saqqaq|Loschbour', files)]
files <- files[grepl('third',files)]
print(files)
dfs <-  do.call(rbind, lapply(files, read.data.frames))
pdf(OUT_path)
print(plotting(dfs))
dev.off()
