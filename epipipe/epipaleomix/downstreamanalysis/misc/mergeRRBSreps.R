mergedataframes <- function(df1,df2){

}

concatdfs <- function(d1, d2){
    mergedataframes(d1,
                   d2)
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

getmean <- function(x,y){
    rowMeans(dfsmerge[,c(x,y)],na.rm=T)
}

files <- list.files(pattern='Bcskeletalmuscle0111002')
dfs <- parallel::mclapply(files,readdf,mc.cores=2)
dfsmerge <- merge(dfs[[1]],dfs[[2]], by='region', suffixes = c(".df1",".df2"),all=T)

finaldf <- data.frame(do.call(rbind,strsplit(dfsmerge$region,'_')), 
                      getmean(11,23),
                      getmean(12,24))
names(finaldf) <- c('chr', 'start', 'end', 'meancoverage', 'meanpercentmethyl')
finaldf$bin <- cut(finaldf$meanpercentmethyl, breaks = c(-1,seq(10, 100,10)),labels=1:10)

f <- function(idx){
mean(finaldf$meanpercentmethy[finaldf$bin==idx],na.rm=T)

}
cbind(1:10,sapply(1:10,f))
## otherstuff
a <- read.table('Saqqaq_MethylMap_methyl450k2000_bedcoord.txt.gz')
b <- a[a$V5>10,]
## Rmin <- min(b$V6[b$V6!=0])
b$V6[b$V6<Rmin] = Rmin
## Rmin <- min(b$V6[b$V6!=0])
## Rmax <- max(b$V6[b$V6!=0])
## jump <- (Rmax-Rmin)/10
## br <-  c(-1,seq(Rmin,Rmax,jump))
br <- c(-1,quantile(b$V6[b$V6!=0],probs=seq(.1,1,.1)))
b$bin <- cut(b$V6, breaks = br,labels=1:10)
f <- function(idx){
mean(b$V6[b$bin==idx],na.rm=T)
}
cbind(1:10,sapply(1:10,f))
