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



round(quantile(vdf$v,probs=seq(0,1,0.1),na.rm=TRUE),2)
###### vdf <- sapply(names(dfs), getvar.gene, SIMPLIFY=TRUE))
with(variancedf, gene[variance<0.02])



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
