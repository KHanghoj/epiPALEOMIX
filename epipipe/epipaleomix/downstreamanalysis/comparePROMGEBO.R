##         data.frame('linearmodel'=n,
##                    'rsquared'=val,
##                    'compared'=paste(mergedf$name.df1[1],mergedf$name.df2[1],sep='_'),
##                    'cov_cpg_cutoff'=paste(cutoffcov, cutoffcpg,sep='_'),
##                    'datapoints'=nrow(mergedf)
##                    )
##     }
## }

## bigfunc <-  function(idx){
##     mergedf <- concatdfs(DFS[[COMBINAT[idx,1]]],
##                           DFS[[COMBINAT[idx,2]]])
##     do.call(rbind, lapply(names(ARGUMENTLIST), runlms, mergedf=mergedf))
## }

readdf <-  function(f, GEBOBOOL){
    df <- read.table(f, colClasses=c('character', 'numeric', 'numeric',
                            'numeric', 'numeric', 'numeric',
                            'numeric', 'numeric', 'numeric'))
    nam <- unlist(strsplit(f, '/'))
    nam <- nam[length(nam)]
    df$name <- unlist(strsplit(nam, '_'))[1]
    if(GEBOBOOL){
        df$GEBOREGION <- with(df, sprintf('%s_%s_%s', V1,V2,V3))
    } else {
        df$PROMREGION <- with(df, sprintf('%s_%s_%s', V1,V2,V3))
    }
    df
}

mergedataframes <- function(idx, dfsG, dfsP){
    prom <- dfsP[[idx]]
    gebo <- dfsG[[idx]]
    mdf <- merge(prom, INFOTABLE, by='PROMREGION')
    mdfall <- merge(gebo, mdf, by='GEBOREGION', suffixes = c(".GEBO",".PROM"))
    ## mdfall <- mdfall[mdfall$PROMCATEGORY=='LOW',]
    ## this is not the optimal way of handling zeros and maximum:::
    mdfall$V6.PROM[mdfall$V6.PROM==0] <- min(mdfall$V6.PROM[mdfall$V6.PROM>0]) 
    maxallowed <- quantile(mdfall$V6.PROM,probs=.95)
    mdfall$V6.PROM[mdfall$V6.PROM>=maxallowed] <- maxallowed
    mdfall$ratio = with(mdfall, V6.GEBO/V6.PROM)
    mdfall$logratio = with(mdfall, log(ratio))
    mdfall$logratio[is.infinite(mdfall$logratio)] =
        with(mdfall, min(logratio[!is.infinite(logratio)]))
    mdfall
}

concatdfs <- function(files){
    GEBOS <- sort(files[grepl('GEBO',files)])
    PROMS <- sort(files[grepl('PROM',files)])
    dfsGEBOS <- lapply(GEBOS,readdf, TRUE)
    dfsPROMS <- lapply(PROMS,readdf, FALSE)
    idx <- 1:length(GEBOS)
    lapply(idx, mergedataframes, dfsG=dfsGEBOS, dfsP=dfsPROMS)
}
fil.cov <- function(cov, cpgs, df){
    df[with(df, V5.GEBO>=cov & V5.PROM>=cov & V8.PROM>=cpgs & V8.GEBO>=cpgs),]
}


bigfunc <- function(idx, dfs){
    dfsone <- dfs[[COMBINAT[idx,1]]]
    dfstwo <- dfs[[COMBINAT[idx,2]]]
    do.call(rbind, lapply(names(ARGUMENTLIST), runlms, dfsone=dfsone,dfstwo=dfstwo))
}

linear.complex <- function(mdf){
    lm(ratio.one~ratio.two*V5.PROM.two*V5.GEBO.two *
       OVERALLPROMCpGCOUNT.two*V8.GEBO.two,   ### *PROMCATEGORY.two*V8.GEBO.two ,
       data=mdf)
}

linear.log.complex <- function(mdf){
    lm(logratio.one~logratio.two*V5.PROM.two*V5.GEBO.two *
       OVERALLPROMCpGCOUNT.two*V8.GEBO.two,   ### *PROMCATEGORY.two*V8.GEBO.two ,
       data=mdf)
}

linear.PROM <- function(mdf){
    lm(V6.PROM.one~V6.PROM.two*OVERALLPROMCpGCOUNT.two*V8.PROM.two, ## PROMCATEGORY.two* ,
       data=mdf)
}
linear.GEBO <- function(mdf){
    lm(V6.GEBO.one~V6.GEBO.two*V8.GEBO.two ,
       data=mdf)
}


runlms <- function(n, dfsone, dfstwo){

    cov <- ARGUMENTLIST[[n]][1]
    cpgs <-  ARGUMENTLIST[[n]][2]
    dfsone <- fil.cov(cov, cpgs, dfsone)
    dfstwo <- fil.cov(cov,cpgs, dfstwo)
    mdf <- merge(dfsone, dfstwo, by='GEBOREGION', suffixes = c(".one",".two"))
    if(nrow(mdf)>100){
        ## l <- linear.PROM(mdf)
        l <- linear.complex(mdf)
        #l <- linear.log.complex(mdf)
        ## l <- linear.GEBO(mdf)
        return(data.frame('linearmodel'='COMPLEX',
                          'rsquared'=summary(l)$r.squared,
                          'compared'=paste(dfsone$name.PROM[1],dfstwo$name.PROM[1],sep='_'),
                          'cov_cpg_cutoff'=paste(cov, cpgs,sep='_'),
                          'datapoints'=nrow(mdf)
                          ))
    }else {return(NULL)}
}



files <- list.files(recursive=T,full.names=T)
## files <- files[grepl('ALL', files)]
files <- files[grepl('bedcoord', files)]
## INFOTABLE = read.table('~/data/bedfiles/PROM_GEBO_autosom_wochr.INFOFILE.FINAL', comment.char='!',h=T)
INFOTABLE = read.table('~/data/bedfiles/PROM_GEBO_HOUSEKEEPING_autosom_wochr.INFOFILE.FINAL',
    comment.char='!',h=T)
testdfs <- concatdfs(files)

COMBINAT <- expand.grid(1:length(testdfs),1:length(testdfs))
COMBINAT <- COMBINAT[COMBINAT[,1]-COMBINAT[,2]!=0,]

ARGUMENTLIST <-  list('none'=c(0,0),
                      '1'=c(10,10),
                      '2'=c(20,20),
                      '3'=c(30,20),
                      '4'=c(50,20),
                      '450'=c(50,20),
                      '5'=c(100,20),
                      '6'=c(150,20),
                      '7'=c(400,20))

mega <- do.call(rbind, parallel::mclapply(1:nrow(COMBINAT),bigfunc, dfs=testdfs,mc.cores=20))
## lapply(1:nrow(COMBINAT),bigfunc, dfs=testdfs)
require(ggplot2)
splot <- function(){
    mega <- do.call(rbind, parallel::mclapply(1:nrow(COMBINAT),bigfunc, dfs=testdfs,mc.cores=20))
    print(ggplot(mega,aes(compared, rsquared, col=cov_cpg_cutoff, size=datapoints))+geom_point() +
    theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
          legend.position='bottom'));dev.off()
}

## require(ggplot2)
## ggplot(testdfs[[1]],aes(PROMCATEGORY,V6.PROM, group=PROMCATEGORY,col=PROMCATEGORY, fill=PROMCATEGORY)) +
##     geom_boxplot();dev.off()
## ggplot(testdfs[[2]],aes(PROMCATEGORY,ratio, group=PROMCATEGORY,col=PROMCATEGORY, fill=PROMCATEGORY)) +
##     geom_boxplot();dev.off()

## write.table(mega, file='crosscomparecutoffboth.txt',row.names=F,col.names=T,quote=F,sep='\t')
## write.table(mega, file='crosscompare.txt',row.names=F,col.names=T,quote=F,sep='\t')


## require(ggplot2)
## a=read.table('crosscompare.txt',h=T)
## a=read.table('crosscomparecutoffboth.txt',h=T)  # this is way more representative, goode.
## b = a[a$datapoints>10000&grepl('complex|semicomplex',a$linearmodel),]
## b = a[a$datapoints>10000&grepl('complex',a$linearmodel),]
## b = a[a$cov_cpg_cutoff=='10_10',]
## b = a[grepl('10_10|20_20', a$cov_cpg_cutoff)&grepl('complex|semicomplex',a$linearmodel),]
## ggplot(b,aes(x=compared, y=rsquared, shape=linearmodel,col=cov_cpg_cutoff, size=datapoints)) +
##     geom_point() + labs(title='minimum 10000 data points')+
##     theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=.5,size=8), legend.position='bottom')

## ggplot(b,aes(x=compared, y=rsquared,col=cov_cpg_cutoff, size=datapoints)) +
##     geom_point() + labs(title='minimum 10000 data points')+
##     theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=.5,size=8), legend.position='bottom')

