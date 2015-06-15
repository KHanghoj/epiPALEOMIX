mergedataframes <- function(df1,df2){

}

concatdfs <- function(d1, d2){
    mergedataframes(d1, d2)
}

readdf <-  function(f){
    df <- read.table(f,skip=1)
    nam <- unlist(strsplit(f, '/'))
    print(table(keep <- df$V5>=10))
    df <- df[keep,]
    df$name <- nam[length(nam)]
    df$region <- with(df, sprintf('%s_%s_%s', V1,V2,V3))
    df
}

getmean <- function(df, x,y){
    rowMeans(df[,c(x,y)],na.rm=T)
}
main <- function(s){
files <- list.files(pattern=s)
dfs <- parallel::mclapply(files,readdf,mc.cores=2)
dfsmerge <- merge(dfs[[1]],dfs[[2]], by='region', suffixes = c(".df1",".df2"),all=T)

finaldf <- data.frame(do.call(rbind,strsplit(dfsmerge$region,'_')), 
                      getmean(dfsmerge, 11,23),
                      getmean(dfsmerge, 12,24))
names(finaldf) <- c('chr', 'start', 'end', 'meancoverage', 'meanpercentmethyl')
finaldf$bin <- cut(finaldf$meanpercentmethyl, breaks = c(-1,seq(10, 100,10)),labels=1:10)
finaldf[,1] = sapply(strsplit(as.character(finaldf[,1]),'chr'),'[[',2)

write.table(finaldf, file=sprintf('%s_comb_wochr.bed',s),row.names=F,col.names=F,quote=F,sep='\t')
}
main('Bcskeletalmuscleh12817n')
main('Osteobl')
main('Bctestisn30')


getmean.three <- function(df, x,y,z){
    rowMeans(df[,c(x,y,z)],na.rm=T)
}


################
#### MELANO THREE replicates
################
files <- list.files(pattern='Melano')
dfs <- parallel::mclapply(files,readdf,mc.cores=3)
dfsmerge <- merge(dfs[[1]],dfs[[2]], by='region', suffixes = c(".df1",".df2"),all=T)
dfsmergefinal <- merge(dfsmerge,dfs[[3]], by='region', suffixes = c(".ding",".dong"),all=T)
finaldf <- data.frame(do.call(rbind,strsplit(dfsmergefinal$region,'_')), 
                      getmean.three(dfsmergefinal, 11,23, 35),
                      getmean.three(dfsmergefinal, 12,24, 36))
names(finaldf) <- c('chr', 'start', 'end', 'meancoverage', 'meanpercentmethyl')
finaldf$bin <- cut(finaldf$meanpercentmethyl, breaks = c(-1,seq(10, 100,10)),labels=1:10)
finaldf[,1] = sapply(strsplit(as.character(finaldf[,1]),'chr'),'[[',2)

write.table(finaldf, file=sprintf('%s_comb_wochr.bed','Melano'),row.names=F,col.names=F,quote=F,sep='\t')

## someterminalstuff to get all bedcoordinates that makes sense
cat *wochr.bed > all.bed
sort -k 1 -k 2n -k 3n all.bed | cut -f 1,2,3|uniq > allcoordinates.bed

library(GenomicRanges)
dat <- read.table('allcoordinates.bed',colClasses=c('character','numeric', 'numeric'))
dat <- dat[dat$V1%in%1:22,]
df <- GRanges(seq=dat$V1, IRanges(start=dat[,2], end=dat[,3]))
keep <-countOverlaps(df, df)==1
dat <- dat[keep,]
write.table(dat, file='allcoordinatestwo.bed',row.names=F,col.names=F,quote=F,sep='\t')

#######
I made a python script to identify cpg sites only. think this is a better approach than the R-script method. Finally, removed cpgs in Y and X chromosome

awk '{print $1"\t"$2-750"\t"$3+750 }' allcoordinatescpgchecked_final.bed > RRBScoordinates.bed
######



################
#### BRAIN only single GOOD
################
files <- list.files(pattern='Bcbrainh11058n')
dfs <- parallel::mclapply(files,readdf,mc.cores=2)
dfsmerge <- sapply(files,readdf)
finaldf <- data.frame(do.call(rbind,strsplit(dfsmerge$region,'_')),dfsmerge$V11, dfsmerge$V12)
names(finaldf) <- c('chr', 'start', 'end', 'meancoverage', 'meanpercentmethyl')
finaldf$bin <- cut(finaldf$meanpercentmethyl, breaks = c(-1,seq(10, 100,10)),labels=1:10)
finaldf[,1] = sapply(strsplit(as.character(finaldf[,1]),'chr'),'[[',2)
write.table(finaldf, file=sprintf('Bcbrainh11058n_comb_wochr.bed',s),row.names=F,col.names=F,quote=F,sep='\t')

## f <- function(idx){
## mean(finaldf$meanpercentmethy[finaldf$bin==idx],na.rm=T)

## }

## cbind(1:10,sapply(1:10,f))


## ## otherstuff
## a <- read.table('Saqqaq_MethylMap_methyl450k2000_bedcoord.txt.gz')
## b <- a[a$V5>10,]
## ## Rmin <- min(b$V6[b$V6!=0])
## b$V6[b$V6<Rmin] = Rmin
## ## Rmin <- min(b$V6[b$V6!=0])
## ## Rmax <- max(b$V6[b$V6!=0])
## ## jump <- (Rmax-Rmin)/10
## ## br <-  c(-1,seq(Rmin,Rmax,jump))
## br <- c(-1,quantile(b$V6[b$V6!=0],probs=seq(.1,1,.1)))
## b$bin <- cut(b$V6, breaks = br,labels=1:10)
## f <- function(idx){
## mean(b$V6[b$bin==idx],na.rm=T)
## }
## cbind(1:10,sapply(1:10,f))
