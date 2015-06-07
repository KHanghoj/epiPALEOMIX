args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
TYPE = args[2]
OUT_path = args[3]
library(data.table)
require(ggplot2)
source('~/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/extendsmooth.R')
#### Rscript ~/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/means.R ~/Desktop/results/allresults/METHYLSAMPLEOUTPUT/ PROM ~/Desktop/

makenames <-  function(f_in, column){
    nam = unlist(strsplit(as.character(f_in),'/'))
    nam = nam[length(nam)]
    as.character(unlist(strsplit(nam,'_'))[column])
}

read.data.frames <-  function(f_in){
    df <- data.frame(fread(sprintf('zcat %s', f_in), skip=1, colClasses=list(character=c(1,5))))
    dfs <-  split(df, df[,5])
    df <- data.frame(do.call(rbind, lapply(dfs, f)))
    colnames(df) = c('val', 'CpGsites', 'totalnumber')
    df$region <- rownames(df)
    df$regionname <- choosename(unlist(strsplit(makenames(f_in, 3), '\\.'))[1])
    df$samplename <- makenames(f_in, 1)
    print(head(df,n=1))
    # print(table(keep <- df$val<quantile(df$val,probs=0.8) & df$totalnumber>=1500 & df$val > 0))
    # print(table(keep <- df$val<quantile(df$val,probs=0.8) & df$totalnumber>=100 & df$val > 0))
    print(table(keep <- df$val<quantile(df$val,probs=0.8) & df$CpGsites>=25))
    df[keep, ]  # keeps the values below the 80 % quantile
}

plotting <- function(df){
    p <- ggplot(df,aes(x=val)) +
        stat_density(aes(ymax = ..density..,  ymin = -..density.., fill = regionname, color = regionname),
                     geom = "ribbon", position = "identity") +
                         facet_grid(samplename ~ regionname)+coord_flip() +
                             theme(legend.position="none", axis.ticks.x = element_blank(),
                                   axis.text.x = element_blank()) +
                                 labs(title=choosename(TYPE), y='', x='Methylation Score')
    p
}

f <- function(df){
    #c(mean(df$V3/df$V4), nrow(df), sum(df$V4))
    c(sum(df$V3)/sum(df$V4), nrow(df), sum(df$V4))
}

files <- list.files(IN_path, pattern='MethylMap',recursive=T,full.names=T)
#files <- files[grepl('third|PROM|GCSKEWED',files)]
files <- files[grepl(TYPE,files)]
#files <- files[grepl('Saqqaq',files)]
print(files)
dfs <-  do.call(rbind, parallel::mclapply(files, read.data.frames, mc.cores=3))

for (n in unique(dfs$samplename)){
    print(n)
    df <- dfs[dfs$samplename==n, ]
    pdf(paste(OUT_path, sprintf('Methylation_Means_%s_%s.pdf',choosename(TYPE), n),sep='/'))
    print(plotting(df))
    dev.off()
}
