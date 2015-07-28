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
    ## df[,col] <- scale(df[,col])  ### might not a good idea to do before PCA. as looking at variance. might just center by mean
    df[,col] <- df[,col]-mean(df[,col])
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

WINDOW <- (c(seq(1,73),74,seq(73,1))) ## weighted filtering window weighing the same in the center
filter.convolution <- function(col, wind=WINDOW, df=mdf){
    filter(col,wind)/sum(wind)
}

files <- list.files(pattern='txt.gz',recursive=T)

files <- files[grepl('_WriteDepth_CONSERVEDARRAY',files)]
files <- files[!grepl('RISE',files)]

dfs <- lapply(files,readdf)
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf, col='depth'),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))
mdf = Reduce(function(x, y) merge(x, y, by=c("nampos"), all=TRUE), lapply(files,rdf, col='score'),accumulate=F)
rcorr(as.matrix(mdf[2:ncol(mdf)]))

mdfnew <- data.frame(pos=mdf[,1], apply(mdf[,2:ncol(mdf)], 2, filter.convolution))

### PCA & HEATMAP
na.to.mean <- function(n){
  M[is.na(M[,n]),n] = colm[n]
  M[,n]
}

M = as.matrix(mdfnew[,-1])
colm <- colMeans(M,na.rm=T)
M.new<-do.call(cbind, parallel::mclapply(names(colm),na.to.mean,mc.cores=4))
colnames(M.new) = names(colm)
M=M.new
rm(M.new)
X<-t(M)%*%M
sites=nrow(M)
X<-X/(sum(diag(X))/(sites-1))
E<-eigen(X)
PC.12 = round(as.numeric(E$values/sum(E$values))[1:2],3)
require(ggplot2)
plotdf = data.frame(E$vectors)
plotdf$nam = colnames(M)
#ggplot(plotdf, aes(X1,X2,col=nam,label=nam))+geom_point()+geom_text()+labs(x=sprintf('PC1 %s %s',PC.12[1]*100, '%'),y=sprintf('PC2 %s %s',PC.12[2]*100,'%'));dev.off()
ggplot(plotdf, aes(X1,X2,col=nam,label=nam))+geom_text()+labs(x=sprintf('PC1 %s %s',PC.12[1]*100, '%'),y=sprintf('PC2 %s %s',PC.12[2]*100,'%'));dev.off()

heatmap(M,  cexCol=0.7, labRow=NA);dev.off()




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



######## TSS analysis of BAM files
df=fread('zcat Saqqaq_WriteDepth_TSS2kb.txt.gz',data.table=F)
dfs = split(df,df[,5]) ## this takes super long
## keep = lapply(dfs,nrow)>1999
ds = df[df$bedcoord=='1_713068_715068',]
ds = cbind(ds, do.call(rbind, strsplit(ds$bedcoord,'_')))

## this part is for regions with less than 1999 bases. so incomplete regions
d=data.frame(genomicpos=713068:715068)
ds = merge(d,ds,by='genomicpos',all=T, fill=T)
ds[is.na(ds)]=0
s = spectrum(ds$score,plot=FALSE)
peaklength = 1/s$freq[which(s$spec==max(s$spec))]
plot(ds$score)
plot(1/s$freq,s$spec, main=peaklength)
dev.off()
##### 

## convolutions of TSS analysis
window <- (c(seq(1,73),74,seq(73,1))) ## weighted filtering window weighing the same in the center
##dsconvolve = convolve(ds$score,window,type='filter')  ## here you remove the first 
dsconvolve = filter(ds$score,window)
## dsconvolve = filter(ds$score,window)/sum(window)  Make the data comparable
dsconvolve[is.na(dsconvolve)] = 0
plot(dsconvolve)
s = spectrum(dsconvolve,plot=FALSE)
peaklength = 1/s$freq[which(s$spec==max(s$spec))]  ## returns the basesize peak
plot(1/s$freq,s$spec, main=peaklength)
dev.off()


### on the read depth
s = spectrum(ds$depth,plot=FALSE)
peaklength = 1/s$freq[which(s$spec==max(s$spec))]
plot(ds$depth)
plot(1/s$freq,s$spec, main=peaklength)
dev.off()

### on the read depth convolved
dsconvolve = filter(ds$depth,window)
dsconvolve[is.na(dsconvolve)] = 0
plot(dsconvolve)
s = spectrum(dsconvolve,plot=FALSE)
peaklength = 1/s$freq[which(s$spec==max(s$spec))]
plot(1/s$freq,s$spec, main=peaklength)
dev.off()




