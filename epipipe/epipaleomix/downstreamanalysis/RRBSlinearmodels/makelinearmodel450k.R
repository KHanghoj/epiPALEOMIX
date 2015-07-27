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
    sample$nampos <- sprintf('%s_%s', sample[,1], sample[,2]+(lengthtocenter))
    output <- data.frame(makemodel(sample, 0, 0),
                         makemodel(sample, 10, 10),
                         makemodel(sample, 20, 10),
                         makemodel(sample, 50, 10),
                         makemodel(sample, 60, 20))

    nam <- unlist(strsplit(f, '/'))
    nam <- nam[length(nam)]
    nam <- unlist(strsplit(nam, '_'))[1]
    modelnames <- c('simplest', 'medium', 'complex')
    totlength = lengthtocenter*2
    colnames(output) <- c(paste(sprintf('%s.%s.0.0',nam,totlength),modelnames, sep='.'),
                          paste(sprintf('%s.%s.10.10',nam,totlength),modelnames, sep='.'),
                          paste(sprintf('%s.%s.20.10',nam,totlength),modelnames, sep='.'),
                          paste(sprintf('%s.%s.50.10',nam,totlength),modelnames, sep='.'),
                          paste(sprintf('%s.%s.60.20',nam,totlength),modelnames, sep='.'))
    output
}
load('complete450kdatasetwithcoord_sorted.Rdata')  #CALLED metdata  samples from column 2 to 57 incl.
metdata$nampos <- with(metdata, sprintf('%s_%s', chr, genomicpos))
metdatanamconverter <- read.table('methyl450kIDtoTissueSample',h=T)
#files <- list.files('temp/', pattern='bedcoord.txt.gz',full.names=T)
files <- list.files('temp/', pattern='Abo',full.names=T)
print(files)
d <- cbind('tissuetype'=metdatanamconverter[,2], do.call(cbind, lapply(files, concatmodels)))
write.table(d, file='lmdata_abo.txt',row.names=F,col.names=T,quote=F,sep='\t')





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

