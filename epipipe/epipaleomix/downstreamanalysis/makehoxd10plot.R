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
    df <- df[df$genomicpos<176998251,]
    ##df <- df[df$genomicpos<176999000,]
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
    
    print(c(sum(df[,3]),sum(df[,4]),nrow(df)))
    df$pos = df$genomicpos
    df$nam <- nam
    #df
    data.frame('pos'=vecpos, 'methylprop'=vecval, 'nam'=nam)
}

getmeans <- function(dfidx){
    df <- DFS[[dfidx]]
    data.frame(pos=df$pos,
               'means'=extendsmooth(df$methylprop,WINSIZElate),  #/
               nam=df$nam[1])
}
WINSIZE = 60
WINHALF = WINSIZE/2
WINSIZElate = 10
files = list.files('HOXD104basesOUTPUT',recursive=T,pattern='.txt.gz',full.names=T)

#files = list.files('HOXD10skipfirstbaseOUTPUT',recursive=T,pattern='.txt.gz',full.names=T)
genes = cbind('nam'=c('hoxd10','hoxd9','hoxd8'),
    data.frame('start'=c(176981492,176987413,176994468),
               'end'=c(176984670,176989645,176997423)))
##chr2:176,979,251-176,998,251
sampleN = sapply(files,getnames)
DFS = lapply(files, readsampledf)
names(DFS) = sampleN
bigdf = do.call(rbind,lapply(names(DFS), getmeans))
require(ggplot2)
p <- ggplot(bigdf,aes(pos,means,col=nam))+geom_line()+
    facet_wrap(~nam, ncol=1, scales='free_y')+
    ##theme(axis.text.x = element_text(angle = 90, hjust = -.5, vjust=0.5, size=8),
            theme(legend.position='None')

pdf('~/Desktop/works.pdf')
## print(ggplot(bigdf,aes(pos,means,col=nam))+geom_line()+facet_wrap(~nam)+
##       theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
##             legend.position='bottom')+
print(p+labs(x='Genomic Position', y='Ms score')+
      geom_segment(aes(x = genes$start[1], y = 0,
                       xend = genes$end[1], yend = 0),col='black')+
      geom_segment(aes(x = genes$start[2], y = 0,
                       xend = genes$end[2], yend = 0),col='black')+
      geom_segment(aes(x = genes$start[3], y = 0,
                       xend = genes$end[3], yend = 0),col='black'));dev.off()
                                        #      facet_wrap(~nam, scales="free_y")+
#dev.off()
