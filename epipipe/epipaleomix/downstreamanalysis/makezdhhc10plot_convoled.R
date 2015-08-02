
getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

readsampledf <- function(f){
    nam <- getnames(f)
    print(nam)
    df <-  read.table(f, comment.char='!',h=T)
    print(table(remove <- (df[,3]/df[,4]>0.33 & df[,3]>3)))
    df <- df[!remove,]

    totdea = sum(df[,3])
    totreads = sum(df[,4])
    CpGsites = nrow(df)
    print(c(totdea, totreads, CpGsites))
    df$pos = df$genomicpos
    newnam <- sprintf("%s;  dea: %s;  totreads: %s;  CpGsites:%s",
                      nam, totdea, totreads, CpGsites)
    methylprop = filter(df[,3]/df[,4], WINDOW)/sum(WINDOW)
    pos = df[!is.na(methylprop),2]
    data.frame('pos'=pos,
               'methylprop'=methylprop[!is.na(methylprop)],
               'nam'=newnam)

}



WINSIZE <- 60
WINDOW <- c(seq(1,(WINSIZE/2)-5),rep(WINSIZE/2, 10), seq((WINSIZE/2)-5,1))
#WINDOW <- rep(1,WINSIZE)
files = list.files(pattern='HOXD10',recursive=T,full.names=T)
print(files)

genes = cbind('nam'=c('ZDHHC4','C7orf26','ZNF853'),
    data.frame('start'=c(6617064, 6629914, 6655526),
               'end'=c(6628610, 6648355, 6663921)))
sampleN = sapply(files,getnames)
bigdf = do.call(rbind, lapply(files, readsampledf))
require(ggplot2)
p <- ggplot(bigdf[bigdf$pos>6606884 & bigdf$pos<6687956 ,],aes(pos,methylprop,col=nam))+geom_line()+
    facet_wrap(~nam, ncol=1, scales='free_y')+
    theme(legend.position='None')+
pdf('ZDHHC10.pdf')
print(p+labs(x='Genomic Position', y='Ms score')+
      geom_segment(aes(x = genes$start[1], y = 0,
                       xend = genes$end[1], yend = 0),col='black')+
      geom_segment(aes(x = genes$start[2], y = 0,
                       xend = genes$end[2], yend = 0),col='black')+
      geom_segment(aes(x = genes$start[3], y = 0,
                       xend = genes$end[3], yend = 0),col='black'));dev.off()
