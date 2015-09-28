f.mostcomplex <- function(mdf){
    ## print('mostcomplex')
    ## print(summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf)))
    summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf))$r.squared
}


## NOT as a covariant Rsquared: 55+: 
f.mostcomplex_no <- function(mdf){
    #print('mostcomplex')
    #print(summary(lm(methylpercent ~ methylprop+cpgread+deaminatedsites+CpGsites+deaminvariance+nodeaminvariance, na.action='na.exclude',data=mdf)))
    summary(lm(methylpercent ~ methylprop+cpgread+deaminatedsites+CpGsites+deaminvariance+nodeaminvariance, na.action='na.exclude',data=mdf))$r.squared
}

getnames <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[1]
}

getnames.unfilt <- function(f){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '\\.'))[1]
}

getlengthtocenter <- function(f){
    nam <- unlist(strsplit(f, '/'))
    as.numeric(unlist(strsplit(unlist(strsplit(nam[length(nam)],'_'))[3],'k'))[2])/2
}

readtissuedf <- function(f){
    df <- fread(f, data.table=FALSE)
    #df <-  read.table(f)
    colnames(df) <- c('chrom', 'pos', 'end', 'cov', 'methylpercent')    
    df$nampos <- sprintf('%s_%s', df[,1], df[,2])
    df$namtis <- getnames.unfilt(f)
    df[,4:7]  ## only take the columns we need later on.
}


readsampledf <- function(f){
    lengthtocenter <- getlengthtocenter(f)
    nam <- getnames(f)
    print(nam)
    df <- fread(sprintf('zcat %s', f), data.table=FALSE,h=T)
    print(nrow(df)) ## removing zeros and top values
    if(nrow(df)>5000){
        df$nampos <- sprintf('%s_%s', df[,1], df[,2]+lengthtocenter)
        df$namsam <- nam
        df
    }
}

makemodel <- function(cov, tissue, sample){
    aux <- which(sample[,5] >= cov) ## it is column 5 after removing GCCONT, as nampos is no longer first column
    sampleselect <- sample[aux,]
    print(c(tissue$namtis[1],sample$namsam[1], cov))
    mdf <- join(sampleselect, tissue, by='nampos',type='inner')
    if(nrow(mdf)>100){data.frame('model'=c('mostcomplex'),
                                 'rsquared'=f.mostcomplex(mdf),
                                 'datapoints'=nrow(mdf),
                                 'cov_cutoff'=cov,
                                 'samplename'=mdf$namsam[1],
                                 'tissuename'=mdf$namtis[1]
                                 )
                  }
}


analpertis <- function(tissuefilespath, sample){
    tissuedf <- readtissuedf(tissuefilespath)
    do.call(rbind, parallel::mclapply(CUTOFFVAL,
                                      makemodel,
                                      tissue=tissuedf,
                                      sample=sample,
                                      mc.cores=12))
}

runanalyses<- function(samplepath){
    sample <- readsampledf(samplepath)
    tissuefiles <- list.files('/home/krishang/data/methylation/RRBS/ucscdownfiltered',pattern='.bed',full.names=T)
#    tissuefiles <- tissuefiles[grepl('Osteo', tissuefiles)]
#    tissuefiles <- tissuefiles[grepl('_1', tissuefiles)]
    do.call(rbind, lapply(tissuefiles, analpertis, sample=sample))
}

require(data.table)
require(plyr)

CUTOFFVAL <- seq(0,250,10)
samplefiles <- list.files('/home/krishang/data/methyl_rawdata/bedcoord_RRBSk1500', pattern='bedcoord.txt.gz',full.names=T)
#samplefiles <- samplefiles[grepl('RRBS', samplefiles)]
#samplefiles <- samplefiles[grepl('Alt', samplefiles)]
print(samplefiles)

mega <- do.call(rbind, lapply(samplefiles, runanalyses))

write.table(mega, file='BIGDF_RRBSlinearmodel_newmethyl_withvariance_ALL.txt',row.names=F,col.names=T,quote=F,sep='\t')

quit('n')


mega <- read.table('RRBSlinearmodel_binomdataeffect.txt',h=T)
mega$compared <- with(mega, paste(samplename, tissuename, sep='_'))
require(ggplot2)
m <- subset(mega, datapoints>5000&!grepl('_2', tissuename)&!grepl('_3', tissuename))
m <- subset(mega,grepl('0_0',cov_cutoff)&grepl('Saqqaq',samplename) &datapoints>50000&grepl('compWCpG',model))
m <- subset(mega,grepl('0_0|10_0',cov_cutoff) &datapoints>50000&grepl('simplest',model))
#m <- subset(mega, cov_cutoff=='0_0' & grepl('simplest|medium',model))
# because is screwed up naming. it is fixed
mod = as.character(mega$model)
mod[seq(2,length(mega$model),5)] = 'simpler'
mega$newmodel = mod
splot <- function(){
      pdf('oneplotnew.pdf')
	print(ggplot(m,aes(compared, rsquared, col=tissuename, shape=cov_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
                	   legend.position='bottom'));dev.off()
}

splot <- function(){
mfirst = m[1:(nrow(m)/2),]
msecond = m[(nrow(m)/2):nrow(m),]
      pdf('oneplotnew.pdf')
	print(ggplot(mfirst,aes(compared, rsquared, col=cov_cutoff, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=5),
                	   legend.position='bottom')+coord_flip()+guides(col=guide_legend(nrow=3)));dev.off()
      pdf('twoplotnew.pdf')
	print(ggplot(msecond,aes(compared, rsquared, col=cov_cutoff, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=5),
                	   legend.position='bottom')+coord_flip()+guides(col=guide_legend(nrow=3)));dev.off()
#guides(col = guide_legend(nrow = 8))
}

m <- subset(mega, datapoints>50000&!grepl('_2', tissuename)&!grepl('_3', tissuename)&grepl('simpler|medium|most',model))

splot <- function(){
midx = seq(nrow(m)/4,nrow(m),nrow(m)/4)
mfirst = m[1:midx[1],]
mmid = m[midx[1]:midx[2],]
mmid2 = m[midx[2]:midx[3],]
msecond = m[midx[3]:midx[4],]
      pdf('1plot.pdf')
	print(ggplot(mfirst,aes(compared, rsquared, col=cov_cutoff, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=5),
                	   legend.position='bottom')+coord_flip()+guides(col=guide_legend(nrow=3)));dev.off()
      pdf('2plot.pdf')
	print(ggplot(mmid,aes(compared, rsquared, col=cov_cutoff, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=5),
                	   legend.position='bottom')+coord_flip()+guides(col=guide_legend(nrow=3)));dev.off()
      pdf('3plot.pdf')
	print(ggplot(mmid2,aes(compared, rsquared, col=cov_cutoff, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=5),
                	   legend.position='bottom')+coord_flip()+guides(col=guide_legend(nrow=3)));dev.off()
      pdf('4plot.pdf')
	print(ggplot(msecond,aes(compared, rsquared, col=cov_cutoff, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=5),
                	   legend.position='bottom')+coord_flip()+guides(col=guide_legend(nrow=3)));dev.off()
#guides(col = guide_legend(nrow = 8))
}

splot <- function(){
      pdf('oneplot.pdf')
	print(ggplot(m,aes(compared, rsquared, col=tissuename, shape=model))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom')+coord_flip());dev.off()
}


splot <- function(){
      pdf('firstpart.pdf')
	print(ggplot(m[1:(nrow(m)/2),],aes(compared, rsquared, col=tissuename, shape=cov_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
      pdf('secondpart.pdf')
	print(ggplot(m[(nrow(m)/2):nrow(m),],aes(compared, rsquared, col=tissuename, shape=cov_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
}

splot <- function(){
      pdf('first.pdf')
	print(ggplot(m[1:300,],aes(compared, rsquared, col=tissuename, shape=cov_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
      pdf('second.pdf')
	print(ggplot(m[300:600,],aes(compared, rsquared, col=tissuename, shape=cov_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()
      pdf('three.pdf')
	print(ggplot(m[600:nrow(m),],aes(compared, rsquared, col=tissuename, shape=cov_cutoff))+geom_point() +
      		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=6),
                	   legend.position='bottom'));dev.off()

}


splot1 <- function(){
    print(ggplot(m,aes(compared, rsquared, col=cov_cutoff, shape=model,size=datapoints))+geom_point() +
          theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
                legend.position='bottom'));dev.off()
}
