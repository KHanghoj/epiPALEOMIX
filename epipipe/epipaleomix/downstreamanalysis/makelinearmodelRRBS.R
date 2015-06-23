f.complex <- function(mdf){
    summary(lm(methylpercent ~ methylprop*coverage, na.action='na.exclude',data=mdf))$r.squared
}
f.simpler <- function(mdf){
    summary(lm(methylpercent ~ methylprop+coverage, na.action='na.exclude',data=mdf))$r.squared
}
f.simplest <- function(mdf){
    summary(lm(methylpercent ~ methylprop, na.action='na.exclude',data=mdf))$r.squared
}
f.compWCpGonly <- function(mdf){
    summary(lm(methylpercent ~ methylprop*coverage*CpG, na.action='na.exclude',data=mdf))$r.squared
}
f.compWCpGacontent <- function(mdf){
    summary(lm(methylpercent ~ methylprop*coverage*GCcont*CpG, na.action='na.exclude',data=mdf))$r.squared
}

getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

getnames.unfilt <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '\\.'))[1]
}

readtissuedf <- function(f){
    ## chrom, pos, end, cov, methylpercent, bin
    ## nam <- getnames(f)
    ## print(nam)
    ## df$nam <- nam
    df <-  read.table(f)
    colnames(df) <- c('chrom', 'pos', 'end', 'cov', 'methylpercent')    
    df$nampos <- sprintf('%s_%s', df[,1], df[,2])
    df$nam <- getnames.unfilt(f)
    print(df$nam[1])
    df
}


readsampledf <- function(f){
    ## chrom, pos, end, cov, methylpercent, bin
    lengthtocenter <- as.numeric(unlist(strsplit(unlist(strsplit(f,'_'))[3],'k'))[2])/2
    nam <- getnames(f)
    print(nam)
    df <-  read.table(f, comment.char='!',h=T)
    print(nrow(df)) ## removing zeros and top values
    # df <- df[df$methylprop>quantile(df$methylprop)[1]&df$methylprop<quantile(df$methylprop)[4],]
    # df <- df[df$methylprop>quantile(df$methylprop)[1],] # remove zeroes only
    if(nrow(df)>5000){
        df$nampos <- sprintf('%s_%s_%s', df[,1], df[,2], df[,3])
        dfm <- merge(df, GCCONT, by='nampos')
        dfm$nampos <- sprintf('%s_%s', dfm[,2], dfm[,3]+lengthtocenter)
        dfm$nam <- nam
        dfm
     }
}

makemodel <- function(argl, tissue, sample){
    cov <- ARGUMENTLIST[[argl]][1]
    cpg <- ARGUMENTLIST[[argl]][2]
    aux <- which(sample[,6] >= cov)
    aux <- aux[which(sample[aux,8] >= cpg)]
    sampleselect <- sample[aux,]
    mdf <- merge(sampleselect,tissue, by='nampos', suffixes=c('sam','tis'))
    if(nrow(mdf)>100){data.frame('model'=c('simplest', 'medium', 'complex','compWCpGonly','compWCpGandcontent'),
                                 'rsquared'=c(f.simplest(mdf), f.simpler(mdf),
                                     f.complex(mdf),f.compWCpGonly(mdf), f.compWCpGacontent(mdf)),
                                 'datapoints'=nrow(mdf),
                                 'cov_cpg_cutoff'=paste(cov, cpg,sep='_'),
                                 'samplename'=mdf$namsam[1],
                                 'tissuename'=mdf$namtis[1]
                                 )
                  }
}

ARGUMENTLIST <-  list('none'=c(0,0),
                      '1a'=c(10,0),
                      '1'=c(10,10),
                      '2'=c(20,20),
                      '3'=c(30,20))
                      ##'4'=c(40,20),
                      ##'5'=c(50,10),
                      ##'6'=c(50,20))
                      ## '3'=c(100,50),
                      ## '4'=c(200,20),
                      ## '450'=c(300,20),
                      ## '5'=c(400,20))

analpertis <- function(tissueidx, sample){
    tissue <- TISSUEDATA[[tissueidx]]
    print(c(tissue$nam[1],sample$nam[1]))
    do.call(rbind, parallel::mclapply(names(ARGUMENTLIST),
                                      makemodel,
                                      tissue=tissue,
                                      sample=sample,
                                      mc.cores=5))
}

concattis <- function(sampleidx){
    sample <- SAMPLEDATA[[sampleidx]]
    ## do.call(rbind, lapply(TISSUE.N, analpertis, sample=sample))
    do.call(rbind, parallel::mclapply(TISSUE.N, analpertis, sample=sample, mc.cores=4))
}

runanalyses <- function(){
    do.call(rbind, lapply(SAMPLE.N, concattis))
}



samplefiles <- list.files('temp', pattern='bedcoord.txt.gz',full.names=T)
samplefiles <- samplefiles[grepl('Altai', samplefiles)]
SAMPLE.N <- sapply(samplefiles, getnames)
print(SAMPLE.N)

GCCONT <- read.table('/home/krishang/data/methylation/RRBS/temp/RRBScoordinates.gccontent')
colnames(GCCONT) <- c('nampos', 'GCcont', 'CpG')
SAMPLEDATA <-  parallel::mclapply(samplefiles, readsampledf, mc.cores=10)
names(SAMPLEDATA) <- SAMPLE.N
SAMPLEDATA <- SAMPLEDATA[!sapply(SAMPLEDATA,is.null)]
SAMPLE.N <- names(SAMPLEDATA)
print(names(SAMPLEDATA))
print(SAMPLE.N)

tissuefiles <- list.files('/home/krishang/data/methylation/RRBS/ucscdownfiltered',pattern='.bed',full.names=T)

TISSUE.N <- sapply(tissuefiles, getnames.unfilt)

TISSUEDATA <-  parallel::mclapply(tissuefiles, readtissuedf, mc.cores=10)
names(TISSUEDATA) <- TISSUE.N

bigdf <- runanalyses()

write.table(bigdf, file='RRBSlinearmodel_binomdataeffect.txt',row.names=F,col.names=T,quote=F,sep='\t')



exit()
mega <- read.table('RRBSlinearmodel_quantilecutoff.txt',h=T)
mega <- read.table('allRRBSlinearmodel.txt',h=T)
mega$compared <- with(mega, paste(samplename, tissuename, sep='_'))
require(ggplot2)
m <- subset(mega, datapoints>5000&grepl('complex|compWCpG',model))
m <- subset(mega,grepl('0_0',cov_cpg_cutoff)&grepl('Saqqaq',samplename) &datapoints>50000&grepl('compWCpG',model))
m <- subset(mega,grepl('0_0|10_0',cov_cpg_cutoff) &datapoints>50000&grepl('simplest',model))
#m <- subset(mega, cov_cpg_cutoff=='0_0' & grepl('simplest|medium',model))

splot <- function(){
	print(ggplot(m,aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
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
