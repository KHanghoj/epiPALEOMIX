f.complex <- function(mdf){
    summary(lm(methylprop.df1 ~ methylprop.df2*coverage.df2, na.action='na.exclude',data=mdf))$r.squared
}
f.simpler <- function(mdf){
    summary(lm(methylprop.df1 ~ methylprop.df2+coverage.df2, na.action='na.exclude',data=mdf))$r.squared
}
f.simplest <- function(mdf){
    summary(lm(methylprop.df1 ~ methylprop.df2, na.action='na.exclude',data=mdf))$r.squared
}
f.compWCpGonly <- function(mdf){
    summary(lm(methylprop.df1 ~ methylprop.df2*coverage.df2*CpG, na.action='na.exclude',data=mdf))$r.squared
}
f.compWCpGacontent <- function(mdf){
    summary(lm(methylprop.df1 ~ methylprop.df2*coverage.df2*GCcont*CpG, na.action='na.exclude',data=mdf))$r.squared
}

getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

readsampledf <- function(f){
    ## chrom, pos, end, cov, methylpercent, cpgcount
    nam <- getnames(f)
    print(nam)
    df <-  read.table(f, comment.char='!',h=T)
    print(nrow(df)) ## removing zeros and top values
    # df <- df[df$methylprop>quantile(df$methylprop)[1]&df$methylprop<quantile(df$methylprop)[4],]
    # df <- df[df$methylprop>quantile(df$methylprop)[1],] # remove zeroes only
    if(nrow(df)>5000){
        df$nampos <- sprintf('%s_%s_%s', df[,1], df[,2], df[,3])
        dfm <- merge(df, GCCONT, by='nampos')
        dfm$nam <- nam
        dfm
     }
}

selectsample <- function(cov, cpg, sample){
    aux <- which(sample[,6] >= cov)
    aux <- aux[which(sample[aux,8] >= cpg)]
    sample[aux,]
}

makemodel <- function(argl, sample1, sample2){
    cov <- ARGUMENTLIST[[argl]][1]
    cpg <- ARGUMENTLIST[[argl]][2]
    sample1 <- selectsample(cov, cpg, sample1)
    sample2 <- selectsample(cov, cpg, sample2)
    mdf <- merge(sample1,sample2, by='nampos', suffixes=c('.df1','.df2'))
    if(nrow(mdf)>100){data.frame('model'=c('simplest', 'medium', 'complex','compWCpGonly','compWCpGandcontent'),
                                 'rsquared'=c(f.simplest(mdf), f.simpler(mdf),
                                     f.complex(mdf),f.compWCpGonly(mdf), f.compWCpGacontent(mdf)),
                                 'datapoints'=nrow(mdf),
                                 'cov_cpg_cutoff'=paste(cov, cpg,sep='_'),
                                 'sample1'=mdf$nam.df1[1],
                                 'sample2'=mdf$nam.df2[1]
                                 )
                  }
}

ARGUMENTLIST <-  list('none'=c(0,0),
                      '1a'=c(10,0),
                      '1'=c(10,5),
                      '2'=c(20,10),
                      '3'=c(30,10))
                      ##'4'=c(40,20),
                      ##'5'=c(50,10),
                      ##'6'=c(50,20))
                      ## '3'=c(100,50),
                      ## '4'=c(200,20),
                      ## '450'=c(300,20),
                      ## '5'=c(400,20))

analpersample <- function(sampleidx){
    sample1 <- SAMPLEDATA[[COMBINANT[sampleidx,1]]]
    sample2 <- SAMPLEDATA[[COMBINANT[sampleidx,2]]]
    print(c(sample1$nam[1],sample2$nam[1]))
    do.call(rbind, parallel::mclapply(names(ARGUMENTLIST),
                                      makemodel,
                                      sample1=sample1,
                                      sample2=sample2,
                                      mc.cores=5))
}

samplefiles <- list.files('chrom1', pattern='bedcoord.txt.gz',full.names=T)
SAMPLE.N <- sapply(samplefiles, getnames)
print(SAMPLE.N)

GCCONT <- read.table('/home/krishang/data/methylation/crosscompare/chrom1/2kchunkschrom1.gccontent')
colnames(GCCONT) <- c('nampos', 'GCcont', 'CpG')

SAMPLEDATA <-  parallel::mclapply(samplefiles, readsampledf, mc.cores=10)
names(SAMPLEDATA) <- SAMPLE.N
SAMPLEDATA <- SAMPLEDATA[!sapply(SAMPLEDATA,is.null)]
SAMPLE.N <- names(SAMPLEDATA)
print(names(SAMPLEDATA))
print(SAMPLE.N)
COMBINANT <- expand.grid(1:length(SAMPLE.N),1:length(SAMPLE.N))
COMBINANT <- COMBINANT[COMBINANT[,1]-COMBINANT[,2]!=0,]
    

bigdf <- do.call(rbind, parallel::mclapply(1:nrow(COMBINANT), analpersample, mc.cores=3))


write.table(bigdf, file='CROSSCOMPARE_CHROM1.txt',row.names=F,col.names=T,quote=F,sep='\t')



exit()

mega <- read.table('RRBSlinearmodel_binomdataeffect.txt',h=T)
mega$compared <- with(mega, paste(samplename, tissuename, sep='_'))
require(ggplot2)
m <- subset(mega, datapoints>5000&grepl('complex|compWCpG',model))
m <- subset(mega,grepl('0_0',cov_cpg_cutoff)&grepl('Saqqaq',samplename) &datapoints>50000&grepl('compWCpG',model))
m <- subset(mega,grepl('0_0|10_0',cov_cpg_cutoff) &datapoints>50000&grepl('simplest',model))
#m <- subset(mega, cov_cpg_cutoff=='0_0' & grepl('simplest|medium',model))

splot <- function(){
      pdf('oneplot.pdf')
	print(ggplot(m,aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
                	   legend.position='bottom'));dev.off()
}

splot <- function(){
      pdf('oneplot.pdf')
	print(ggplot(m,aes(compared, rsquared, col=tissuename, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
                	   legend.position='bottom'));dev.off()
}


splot <- function(){
      pdf('firstpart.pdf')
	print(ggplot(m[1:(nrow(m)/2),],aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
      pdf('secondpart.pdf')
	print(ggplot(m[(nrow(m)/2):nrow(m),],aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
}

splot <- function(){
      pdf('first.pdf')
	print(ggplot(m[1:300,],aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
      pdf('second.pdf')
	print(ggplot(m[300:600,],aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
      pdf('three.pdf')
	print(ggplot(m[600:nrow(m),],aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()

}


splot1 <- function(){
    print(ggplot(m,aes(compared, rsquared, col=cov_cpg_cutoff, shape=model,size=datapoints))+geom_point() +
          theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
                legend.position='bottom'));dev.off()
}
