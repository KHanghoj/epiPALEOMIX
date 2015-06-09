## saq <- read.table('Saqqaq_MethylMap_PROMHIGH_bedcoord.txt.gz')
## ust <- read.table('UstIshim_MethylMap_PROMHIGH_bedcoord.txt.gz')
## saq$V10 <- with(saq, sprintf('%s_%s_%s', V1,V2,V3))
## ust$V10 <- with(ust, sprintf('%s_%s_%s', V1,V2,V3))
## df <- merge(saq,ust, by='V10', suffixes = c(".saq",".ust"))
## df$delta <- with(df,V5.saq-V5.ust)

## alt <- read.table('/Users/kehanghoej/Desktop/AltaiNeanderthal_MethylMap_PROMHIGH_bedcoord.txt.gz')
## den <- read.table('/Users/kehanghoej/Desktop/Denisova_MethylMap_PROMHIGH_bedcoord.txt.gz')
## alt$V10 <- with(alt, sprintf('%s_%s_%s', V1,V2,V3))
## den$V10 <- with(den, sprintf('%s_%s_%s', V1,V2,V3))
## df <- merge(alt,den, by='V10', suffixes = c(".alt",".den"))
## df$delta <- with(df,V5.alt-V5.den)

## summary(lm(V9.alt~V9.den*delta  , data=df))

f.complex <- function(x, bigdf, df){
    summary(lm(bigdf[,x] ~ df$V6*df$V5, na.action='na.exclude'))$r.squared
}
f.simpler <- function(x, bigdf, df){
    summary(lm(bigdf[,x] ~ df$V6+df$V5, na.action='na.exclude'))$r.squared
}
f.simplest <- function(x, bigdf, df){
    summary(lm(bigdf[,x] ~ df$V6, na.action='na.exclude'))$r.squared
}

makemodel <- function(sample, cutoffcov, cutoffcpg){

    aux <- which(sample[,'V5'] >= cutoffcov)
    aux <- aux[which(sample[aux,'V8'] >= cutoffcpg)]
    sampleselect <- sample[aux,]

    print(table(keep <- metdata$nampos %in% sampleselect$nampos))
    bigdf <- metdata[keep,2:57]
    print(sprintf('%s_%s',nrow(bigdf), nrow(sampleselect)))
    cbind(do.call(rbind, parallel::mclapply(1:ncol(bigdf),f.simplest,
                                            bigdf=bigdf, df=sampleselect, mc.cores=25)),
          do.call(rbind, parallel::mclapply(1:ncol(bigdf),f.simpler,
                                            bigdf=bigdf, df=sampleselect, mc.cores=25)),
          do.call(rbind, parallel::mclapply(1:ncol(bigdf),f.complex,
                                            bigdf=bigdf, df=sampleselect, mc.cores=25)))
}

concatmodels <- function(f){
    lengthtocenter <- as.numeric(unlist(strsplit(unlist(
        strsplit(f,'_'))[3],'k'))[2])/2
    print(f)
    
    sample <- read.table(f)
    sample$nampos <- sprintf('%s_%s', sample[,1], sample[,2]+lengthtocenter)
    ## output <- data.frame(metdatanamconverter[,2],
    ##                      makemodel(sample, 100, 10),
    ##                      makemodel(sample, 100, 50),
    ##                      makemodel(sample, 200, 50),
    ##                      makemodel(sample, 300, 50),
    ##                      makemodel(sample, 400, 50))
    output <- data.frame(makemodel(sample, 100, 10),
                         makemodel(sample, 100, 50),
                         makemodel(sample, 200, 50),
                         makemodel(sample, 300, 50),
                         makemodel(sample, 400, 50))

    nam <- unlist(strsplit(f, '/'))
    nam <- nam[length(nam)]
    nam <- unlist(strsplit(nam, '_'))[1]
    modelnames <- c('simplest', 'medium', 'complex')
    colnames(output) <- c(paste(sprintf('%s.100.10',nam),modelnames, sep='_'),
                          paste(sprintf('%s.100.50',nam),modelnames, sep='_'),
                          paste(sprintf('%s.200.50',nam),modelnames, sep='_'),
                          paste(sprintf('%s.300.50',nam),modelnames, sep='_'),
                          paste(sprintf('%s.400.50',nam),modelnames, sep='_'))
    output
}
#f='temp/Saqqaq_MethylMap_ALLPROM_hashremoved_bedcoord.txt'
#f <- 'temp/Saqqaq4bases_MethylMap_methyl450k1500_bedcoord.txt.gz'
load('complete450kdatasetwithcoord_sorted.Rdata')  #CALLED metdata  samples from column 2 to 57 incl.
metdata$nampos <- with(metdata, sprintf('%s_%s', chr, genomicpos))
metdatanamconverter <- read.table('methyl450kIDtoTissueSample',h=T)
files <- list.files('temp/', pattern='bedcoord.txt.gz',full.names=T)
print(files)
d <- cbind(metdatanamconverter[,2], do.call(cbind, lapply(files, concatmodels)))
write.table(d, file='lmdata.txt',row.names=F,col.names=T,quote=F,sep='\t')





## f.max <- function(x){
##     max(bigdf[,x],na.rm=T)
## }
## sapply(1:ncol(bigdf),f.max)



## Getrows <- function(rowidx, Mat, ancsample){
##     row <- unlist(ancsample[rowidx,1:3])
##     aux <- which(Mat[,'chr'] == row[1])
##     aux <- aux[which(Mat[aux,'genomicpos'] >= row[2])]
##     aux <- aux[which(Mat[aux,'genomicpos'] < row[3])]
##     if(length(aux)>0){c(row, length(aux),colMeans(Mat[aux, 2:57]))}
## }

## dfs <- do.call(rbind, parallel::mclapply(1:nrow(dfsaq), getrows,
##                                          Mat=metdata,
##                                          ancsample=dfsaq, mc.cores=3))
## rownames(dfs) <- sprintf('%s_%s_%s', dfs[,1], dfs[,2], dfs[,3])

## dfsaq$nam <-sprintf('%s_%s_%s', dfsaq[,1], dfsaq[,2], dfsaq[,3])
## dfsaq <- dfsaq[!duplicated(dfsaq$nam),]

## keep <- dfsaq$nam %in% rownames(dfs)
## dfsaq <- dfsaq[keep,]



## dfsnew <- sapply(dfs[, 6:ncol(dfs)], as.numeric)
## dim(dfs)
## dim(a)


## cor(cbind(dfs,a$V9))
## 1  860029  861529

