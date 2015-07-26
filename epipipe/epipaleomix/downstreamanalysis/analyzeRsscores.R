library(data.table)
df=fread('HKjoinedwithRs_noAbo.txt',data.table=FALSE)
dfs <- split(df,df[,1])
names(dfs) <- sort(unique(df[,1]))

getvar.gene <- function(nam){
  d <- dfs[[nam]]
  if(nrow(d)>4){
    data.frame(n=nam,v=var(d$V28),m=mean(d$V28), no=nrow(d))
  }
  
}
vdf <- do.call(rbind, (lapply(names(dfs), getvar.gene)))
vdf <- vdf[order(vdf$n),]  ## order by genename

write.table(vdf$n[quantile(vdf$v,probs=0.9)<=vdf$v], file='HKgene_topvar.txt',row.names=F,col.names=F,quote=F,sep='\t')
write.table(vdf$n[quantile(vdf$v,probs=0.1)>=vdf$v], file='HKgene_lowvar.txt',row.names=F,col.names=F,quote=F,sep='\t')
## correlation studies between samples:
library(gtools)
IDXES <- combinations(length(names(dfs)),2)

### df <- df[,c(1,24,26)]
##### CORRELATION STUDIES
## this is housekeeping
df <- fread('HKjoinedwithRs_withoutABO.txt',data.table=FALSE)
dfs <- split(df[,c(15,28)],df[,14])  ## Normalized data # normalized is better
names(dfs) <- sort(unique(df[,14]))

## this is all genes
df <- fread('zcat Rsscores2000_withoutABO.txt.gz',data.table=F)
## dfs <- split(df[,c(13,26)],df[,12])  ## raw data
dfs <- split(df[,c(13,27)],df[,12])  ## Normalized data # normalized is better
## dfs <- split(df[,c(13,19)],df[,12])  ## methylPROPprom ONLY
names(dfs) <- sort(unique(df[,12]))
renewnames <- function(nam){
  d <- dfs[[nam]]
  colnames(d) <- c('namposPROM', nam)
  d
}
require(Hmisc)
dfs <- lapply(names(dfs),renewnames)
data <- Reduce(function(x, y) merge(x, y, by=c("namposPROM"),all=TRUE), dfs, accumulate=F)
rcorr(as.matrix(data[2:ncol(data)]))


### PCA & HEATMAP
na.to.mean <- function(n){
  M[is.na(M[,n]),n] = colm[n]
  M[,n]
}

M = as.matrix(data[,-1])
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

## canonical correlation analyses
## Altai+Denisova vs REST
## Altai+Denisova vs REST (-Motala -Kostenki) as somewhat less data there
## Altai+Denisova vs REST (-Motala ) as little data there
## Kostenki+UstIshim vs Loschbour|Stuttgart (Late Pleistocene vs Holocene)
## Kostenki+UstIshim vs Loschbour vs Stuttgart (Late Pleistocene vs preNeolithic vs early farmer)
## Kostenki+UstIshim+Loschbour (i.e preNeolithic) vs Stuttgart (early farmer)
## Saqqaq|Aborigine vs REST (ie looking for tissue differences)
## [1] "AltaiNeanderthal" "Denisova"         "Kostenki"         "Loschbour"
## [5] "Motala"           "Saqqaq"           "Stuttgart"        "UstIshim"
require(CCA)
canon.correlation <- function(l, df=M){
  df.left <- df[,l$left]
  df.right <- df[,l$right]
  # display the canonical correlations in out1
  out1 <- cc(df.left,df.right)
  # compute canonical loadings
  list('cccor'=out1, 'ccload'=comput(df.left, df.right, out1))
}
canon.correlation <- function(l, df=M){
  df.left <- df[,l$left]
  df.right <- df[,l$right]
  # display the canonical correlations in out1
  cc(df.left,df.right)
}


canon.matrix.correlation <- function(l, df=M){
  df.left <- df[,l$left]
  df.right <- df[,l$right]
  matcor(df.left,df.right) 
}


canon.correlation(list('left'=c('AltaiNeanderthal', 'Denisova'), 'right'=c('Kostenki','Loschbour', 'Motala', 'Saqqaq', 'Stuttgart', 'UstIshim')))
canon.correlation(list('left'=c('AltaiNeanderthal', 'Denisova'), 'right'=c('Kostenki','Loschbour', 'Saqqaq', 'Stuttgart', 'UstIshim')))
canon.correlation(list('left'=c('AltaiNeanderthal', 'Denisova'), 'right'=c('Loschbour', 'Saqqaq', 'Stuttgart', 'UstIshim')))
canon.correlation(list('left'=c('Kostenki','UstIshim'), 'right'=c('Loschbour')))
canon.correlation(list('left'=c('Kostenki','UstIshim'), 'right'=c('Stuttgart')))
canon.correlation(list('left'=c('Kostenki','UstIshim'), 'right'=c('Stuttgart','Loschbour')))
canon.correlation(list('left'=c('Kostenki','UstIshim','Loschbour'), 'right'=c('Stuttgart')))
## list('left'=c('Kostenki','UstIshim'), 'right'=c('Stuttgart','Loschbour')) can i compare three???
canon.correlation(list('left'=c('Saqqaq'), 'right'=c('AltaiNeanderthal', 'Denisova', 'Kostenki','Loschbour', 'Motala', 'Stuttgart', 'UstIshim')))

M = as.matrix(data[,-1])
colm <- colMeans(M,na.rm=T)
M.new<-do.call(cbind, parallel::mclapply(names(colm),na.to.mean,mc.cores=4))
colnames(M.new) = names(colm)
M=M.new
rm(M.new)



vM = apply(M,1,var)
df[(df$namposPROM %in% data$namposPROM[quantile(vM,probs=0.99)<vM]),'ID']  # the genes showing greatest variance in Rs
write.table(unique(df[(df$namposPROM %in% data$namposPROM[quantile(vM,probs=0.99)<vM]),'ID']),file='greatestvariance_geneID.txt',row.names=F,col.names=F,quote=F,sep='\t')

round(quantile(vdf$v,probs=seq(0,1,0.1),na.rm=TRUE),2)
###### vdf <- sapply(names(dfs), getvar.gene, SIMPLIFY=TRUE))
with(variancedf, gene[variance<0.02])


library(data.table)
df <- fread('Rsscorefortest_withoutHK.txt',data.table=FALSE)
dfs <- split(df,df[,12])
names(dfs) <- sort(unique(df[,12]))

calc.top.bottom <- function(nam){
  d <- dfs[[nam]]
  out <- d$ID[quantile(d$RsSCORES,probs=0.90)-1e-7<=d$RsSCORES]
  write.table(out[order(out)], file=sprintf('gene_%s_top.txt',nam),row.names=F,col.names=F,quote=F,sep='\t')
  out <- d$ID[quantile(d$RsSCORES,probs=0.10)+1e-7>=d$RsSCORES]
  write.table(out[order(out)], file=sprintf('gene_%s_lower.txt',nam),row.names=F,col.names=F,quote=F,sep='\t')
}

sapply(names(dfs), calc.top.bottom)
