require(data.table)
require(Hmisc)
readdf <- function(f){
    print(getnames(f))
    df <- fread(spriBntf('zcat %s',f),data.table=F)
    df$nampos <- sprintf('%s_%s',df[,1],df[,2])
    df$nam <- getnames(f)
    df
}


rdf <- function(f, col){
    print(getnames(f))
    df <- fread(sprintf('zcat %s',f),data.table=F)
    df$nampos <- sprintf('%s_%s',df[,1],df[,2])
    df <- df[,c(col,'nampos')]
    colnames(df) <- c(getnames(f), 'nampos')
    df
}
spectrumwrapper <- function(score){
    ## takes a vector of scores as input
    d <- spectrum(scores,method='pgram',plot=F)
    data.frame(x=1/d$freq,
               y=d$spec)
}
###
###OR
###
####  spectrum(na.omit(mdf)[,2:ncol(mdf)],plot=F)  ## this works too

getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

files <- list.files(pattern='txt.gz',recursive=T)

files <- files[grepl('_WriteDepth_CONSERVEDARRAY',files)]
files <- files[!grepl('RISE',files)]
dfs <- lapply(files,readdf)
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf, col='depth'),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf, col='score'),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))





## correlation and convolution:
## http://www.cs.umd.edu/~djacobs/CMSC426/Convolution.pdf
## https://en.wikipedia.org/wiki/Convolution


mdfnaomit = naomit(mdf)

filteredsample <- function(nam, mdf=mdfnaomit, window=window){
    centered <- (mdfnaomit[,nam]-mean(mdfnaomit[,nam]))/sd(mdfnaomit[,nam])
    ###filter(centered,window)
    convolve(centered,window,type='filter')
}

centeredsaq <- (mdfnaomit[,c('Saqqaq')]-mean(mdfnaomit[,c('Saqqaq')]))/sd(mdfnaomit[,c('Saqqaq')])
centeredabo <- (mdfnaomit[,c('Aborigin')]-mean(mdfnaomit[,c('Aborigin')]))/sd(mdfnaomit[,c('Aborigin')])


window <- (c(seq(1,73),74,seq(73,1))) ## weighted filtering window weighing the same in the center
windowones=rep(1,147)  ## this is not good. just for testing

### Filter:
valsaq <- convolve(centeredsaq,window,type='filter')
valabo <- convolve(centeredabo,window,type='filter')
### OR:
valsaq <- filter(centeredsaq,window)  ## filter is the same as convolve(...,type='filter'). though it adds NA' (window/2) as beginning and end.
valabo <- filter(centeredabo,window)

valdeni <- filteredsample('Denisova')

#co1 <- convolve(datsaq,databo)

jump=5000
for (idx in seq(0,length(valsaq)-jump,jump)){
    dat = ccf(valabo[idx:(idx+jump)],valsaq[idx:(idx+jump)],plot=F)
    l <- which(dat$acf==max(dat$acf))
    print(c(round(max(dat$acf),1),dat$lag[l],idx, idx+jump))
}
    
err.func <- function(p,X1, X2){
    sum(sqrt((X1-X2+p[1])**2))
}
err.func <- function(X1, X2){
    sum(sqrt((X1-X2)**2))/length(X1)
}
### this does not work
## err.func <- function(p,X1, X2){
##     cor(X1,X2+p[1])
## }

dat = ccf(valabo[1:10000],valsaq[1:10000],plot=F)
dat = ccf(valabo,valsaq,plot=F)
par <- c(0,0)
optimx(par,err.func, X1=valsaq,X2=valabo)
##optimx(par,err.func, X1=mdfnaomit[1:50000,c('Saqqaq')],X2=mdfnaomit[1:50000,c('Aborigin')])


cor.test(mdfnaomit[,c('Aborigin')], mdfnaomit[,c('Saqqaq')])
cor.test(valsaq,valabo)


plot(valabo[1:2000])
plot(mdfnaomit[73:2073,c('Aborigin')])
dev.off()

convolved = convolve(mdfnaomit[1:2000,c('Aborigin')],mdfnaomit[1:2000,c('Aborigin')],type='open')
which(convolved==max(convolved))  # gives 2000 ## nothing to skew as same sample.
convolved = convolve(mdfnaomit[1:500,c('Aborigin')],mdfnaomit[1:500,c('Saqqaq')],type='open')
convolved = convolve(valabo[1:500],valsaq[1:500],type='open')
which(convolved==max(convolved))  # gives 492. skew 8 bases

convolved = convolve(mdfnaomit[1:500,c('Aborigin')],mdfnaomit[1:500,c('Saqqaq')],type='open')
convolved = convolve(valabo[1:500],valsaq[1:500],type='open')
which(convolved==max(convolved))  # gives 492. skew 8 bases

co1 = convolve(mdfnaomit[1:300,c('Aborigin')],mdfnaomit[1:300,c('Saqqaq')],type='open')
plot(co1)
which(co1==max(co1))
co2 = convolve(valabo[1:300],valsaq[1:300],type='open')
plot(co2)
which(co2==max(co2))
dev.off()



valsaq = convolve(mdfnaomit[,c('Saqqaq')],window,type='filter')
valabo = convolve(mdfnaomit[,c('Denisova')],window,type='filter')
plot(valsaq[1:500])
plot(valabo[1:500])
dev.off()
