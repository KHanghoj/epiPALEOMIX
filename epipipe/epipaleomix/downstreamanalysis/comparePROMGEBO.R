mergedataframes <- function(df1,df2){
    mdf <- merge(df1,df2, by='V10', suffixes = c(".df1",".df2"))
    print(nrow(mdf))
    merge(mdf, gccont, by.x='V10', by.y='region')
}
concatdfs <- function(d1, d2){
    mergedataframes(d1,
                   d2)
}

runlms <- function(idx, mergedf){
    cutoffcov = ARGUMENTLIST[[idx]][1]
    cutoffcpg = ARGUMENTLIST[[idx]][2]

                                        # mergedf <-  mergedf[mergedf$V5.df1>=cutoffcov & mergedf$V8.df1>=cutoffcpg,]
    mergedf <-  mergedf[mergedf$V5.df1>=cutoffcov & mergedf$V8.df1>=cutoffcpg &
                        mergedf$V5.df2>=cutoffcov & mergedf$V8.df2>=cutoffcpg ,]
    n <- c('simplest', 'simpler', 'semicomplex', 'complex')
    if(nrow(mergedf)>10){
        val <- c(summary(f.simplest(mergedf))$r.squared,
                 summary(f.simpler(mergedf))$r.squared,
                 summary(f.complex.1(mergedf))$r.squared,
                 summary(f.complex(mergedf))$r.squared
                 )        
        data.frame('linearmodel'=n,
                   'rsquared'=val,
                   'compared'=paste(mergedf$name.df1[1],mergedf$name.df2[1],sep='_'),
                   'cov_cpg_cutoff'=paste(cutoffcov, cutoffcpg,sep='_'),
                   'datapoints'=nrow(mergedf)
                   )
    }
}

bigfunc <-  function(idx){
    mergedf <- concatdfs(DFS[[COMBINAT[idx,1]]],
                          DFS[[COMBINAT[idx,2]]])
    do.call(rbind, lapply(names(ARGUMENTLIST), runlms, mergedf=mergedf))
}

readdf <-  function(f){
    df <- read.table(f)
    nam <- unlist(strsplit(f, '/'))
    nam <- nam[length(nam)]
    df$name <- unlist(strsplit(nam, '_'))[1]
    df$V10 <- with(df, sprintf('%s_%s_%s', V1,V2,V3))
    df
}

gccont <- read.table('methyl450k_2000_wochr.gccontentnew', h=T)
files <- list.files(recursive=T,full.names=T)

grepl('PROM|GEBO', files)

print(files)

library(gtools)
#COMBINAT <- combinations(length(files),2, set=TRUE, repeats.allowed=FALSE)
COMBINAT <- expand.grid(1:length(files),1:length(files))
COMBINAT <- COMBINAT[COMBINAT[,1]-COMBINAT[,2]!=0,]
ARGUMENTLIST <-  list('none'=c(0,0),
                      '1'=c(10,10),
                      '2'=c(20,20),
                      '3'=c(30,20),
                      '4'=c(50,20),
                      '450'=c(50,50),
                      '5'=c(100,20),
                      '6'=c(150,30))

DFS <- parallel::mclapply(files,readdf, mc.cores=10)
names(DFS) = 1:length(files)
mega <- do.call(rbind, parallel::mclapply(1:nrow(COMBINAT),bigfunc, mc.cores=20))
write.table(mega, file='crosscomparecutoffboth.txt',row.names=F,col.names=T,quote=F,sep='\t')
write.table(mega, file='crosscompare.txt',row.names=F,col.names=T,quote=F,sep='\t')


require(ggplot2)
a=read.table('crosscompare.txt',h=T)
a=read.table('crosscomparecutoffboth.txt',h=T)  # this is way more representative, goode.
b = a[a$datapoints>10000&grepl('complex|semicomplex',a$linearmodel),]
b = a[a$datapoints>10000&grepl('complex',a$linearmodel),]
b = a[a$cov_cpg_cutoff=='10_10',]
b = a[grepl('10_10|20_20', a$cov_cpg_cutoff)&grepl('complex|semicomplex',a$linearmodel),]
ggplot(b,aes(x=compared, y=rsquared, shape=linearmodel,col=cov_cpg_cutoff, size=datapoints)) +
    geom_point() + labs(title='minimum 10000 data points')+
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=.5,size=8), legend.position='bottom')

ggplot(b,aes(x=compared, y=rsquared,col=cov_cpg_cutoff, size=datapoints)) +
    geom_point() + labs(title='minimum 10000 data points')+
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=.5,size=8), legend.position='bottom')

