require(data.table)
require(Hmisc)
rdf <- function(f, col='deaminatedsites'){
    df <- fread(sprintf('zcat %s',f),data.table=F)
    df$nampos <- sprintf('%s_%s_%s',df[,1],df[,2],df[,3])
    ## df$meth <- df$deaminatedsites/df$coverage
    ## df <- df[,c('meth','nampos')]
    ## df <- df[,c('methylprop','nampos')]
    df <- df[,c(col,'nampos')]
    ## df <- df[,c('coverage','nampos')]
    colnames(df) <- c(getnames(f), 'nampos')
    df
}
getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

files <- list.files(pattern='txt.gz')
#files <- files[grepl('_PROM_',files)]
files <- files[grepl('_GEBO_',files)]
#altai <- rdf('AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz')
#deni <- rdf('Denisova_MethylMap_PROM_bedcoord.txt.gz')
#saq <- rdf('Saqqaq_MethylMap_PROM_bedcoord.txt.gz')

mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))

