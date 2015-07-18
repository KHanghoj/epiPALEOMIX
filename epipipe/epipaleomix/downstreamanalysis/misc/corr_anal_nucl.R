require(data.table)
require(Hmisc)
rdf <- function(f, col){
    print(getnames(f))
    df <- fread(sprintf('zcat %s',f),data.table=F)
    df$nampos <- sprintf('%s_%s_%s',df[,1],df[,2],df[,3])
    df <- df[,c(col,'nampos')]
    colnames(df) <- c(getnames(f), 'nampos')
    df
}
getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

files <- list.files(pattern='txt.gz',recursive=T)

files <- files[grepl('_WriteDepth_CONSERVEDARRAY',files)]

mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf, col='depth'),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf, col='score'),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))

