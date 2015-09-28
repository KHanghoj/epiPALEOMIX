##Rscript  /Users/kehanghoej/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/makehoxd10plot.R

smoothfunc <- function(idx, seq, windowsize){
    sum(seq[idx:(idx+windowsize-1)])/(windowsize)
}

extendsmooth <- function(seq, windowsize=2){
    idx = 1:length(seq)
    smooth = sapply(idx, smoothfunc, seq=seq, windowsize=windowsize)
    keep = !is.na(smooth)
    smooth = smooth[keep]
    additions = sum(!keep)
    if (additions %% 2 == 0){
        smooth = c(rep(smooth[1],additions/2),smooth,
            rep(smooth[length(smooth)],additions/2))
    } else {
        smooth = c(smooth, smooth[length(smooth)])
        additions = additions-1
        smooth = c(rep(smooth[1],additions/2),smooth,
            rep(smooth[length(smooth)],additions/2))
    }
    smooth
}

getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}
readsampledf <- function(f){
    nam <- getnames(f)
    print(nam)
    df <-  read.table(f, comment.char='!',h=T)
    vecpos = c()
    vecval = c()
    for (idx in seq(WINHALF,(nrow(df)-WINHALF),1)){
        start <- idx-WINHALF
        end <- idx+WINHALF
        vecval <- c(vecval,sum(df[start:end,3],na.rm=T)/sum(df[start:end,4],na.rm=T))
        ## vecval <- c(vecval,mean(df[start:end,3]/df[start:end,4],na.rm=T))
        vecpos <- c(vecpos, mean(df[start:end,2]))
    }
    df$methylprop = df[,3]/df[,4]
    totdea = sum(df[,3])
    totreads = sum(df[,4])
    CpGsites = nrow(df)
    print(c(totdea, totreads, CpGsites))
    df$pos = df$genomicpos
    newnam <- sprintf("%s;  dea: %s;  totreads: %s;  CpGsites:%s",
                      nam, totdea, totreads, CpGsites)

    data.frame('pos'=vecpos,
               'methylprop'=vecval,
               'nam'=newnam)
}

getmeans <- function(dfidx){
    df <- DFS[[dfidx]]
    ## data.frame(pos=df$pos,
    ##            'means'=extendsmooth(df$methylprop,WINSIZElate),  #/
    ##            nam=df$nam[1]
    df$means=extendsmooth(df$methylprop,WINSIZElate)
    df
}
WINSIZE = 70
WINHALF = WINSIZE/2
WINSIZElate = 1
#files = list.files('HOXD104basesOUTPUT',recursive=T,pattern='.txt.gz',full.names=T)
files = list.files(pattern='conserved',recursive=T,full.names=T)
files = files[grepl('.txt.gz$', files)]
print(files)

sampleN = sapply(files,getnames)
DFS = lapply(files, readsampledf)
names(DFS) = sampleN
bigdf = do.call(rbind,lapply(names(DFS), getmeans))
require(ggplot2)
p <- ggplot(bigdf,aes(pos,means,col=nam))+geom_line()+
    facet_wrap(~nam, ncol=1, scales='free_y')+
    theme(legend.position='None')

pdf('Conservedregion.pdf')
print(p+labs(x='Genomic Position', y='Ms score'))
dev.off()
