# f.complex <- function(mdf){
#     summary(lm(methylpercent ~ cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf))$r.squared
# }
# f.simpler <- function(mdf){
#     summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites, na.action='na.exclude',data=mdf))$r.squared
# }
# f.simplest <- function(mdf){
#     summary(lm(methylpercent ~ methylprop*cpgread, na.action='na.exclude',data=mdf))$r.squared
# }


## f.mediumcov <- function(mdf){
##     summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance, na.action='na.exclude',data=mdf))$r.squared
## }
## f.mostcomplexcov <- function(mdf){
##     summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf))$r.squared
## }
f.medium <- function(mdf){
    # print('medium')
    # print(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance, na.action='na.exclude',data=mdf))
    summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance, na.action='na.exclude',data=mdf))$r.squared
}

## as a covariant ost Rsquared: 60+:
f.mostcomplex <- function(mdf){
    print('mostcomplex')
    print(summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf)))
    summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf))$r.squared
}
## NOT as a covariant Rsquared: 55+: 
f.mostcomplex <- function(mdf){
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

readtissuedf <- function(f){
    df <-  read.table(f)
    colnames(df) <- c('chrom', 'pos', 'end', 'cov', 'methylpercent')    
    df$nampos <- sprintf('%s_%s', df[,1], df[,2])
    df$nam <- getnames.unfilt(f)
    #print(df$nam[1])
    df
}


readsampledf <- function(f){
    lengthtocenter <- as.numeric(unlist(strsplit(unlist(strsplit(f,'_'))[3],'k'))[2])/2
    nam <- getnames(f)
    print(nam)
    df <-  read.table(f, comment.char='!',h=T)
    print(nrow(df)) ## removing zeros and top values
    if(nrow(df)>5000){
        df$nampos <- sprintf('%s_%s', df[,1], df[,2]+lengthtocenter)
        df$nam <- nam
        df
     }
}

makemodel <- function(cov, tissue, sample){
    aux <- which(sample[,5] >= cov) ## it is column 5 after removing GCCONT, as nampos is no longer first column
    ## aux <- which(sample[,6] >= cov)
    sampleselect <- sample[aux,]
    print(c(tissue$nam[1],sample$nam[1], cov))
    mdf <- merge(sampleselect,tissue, by='nampos', suffixes=c('sam','tis'))
    if(nrow(mdf)>100){data.frame('model'=c('medium','mostcomplex'),
                                 'rsquared'=c(f.medium(mdf), f.mostcomplex(mdf)),
                                 'datapoints'=nrow(mdf),
                                 'cov_cutoff'=cov,
                                 'samplename'=mdf$namsam[1],
                                 'tissuename'=mdf$namtis[1]
                                 )
                  }
}
CUTOFFVAL <- seq(0,250,10)
#CUTOFFVAL <- seq(0,250,20)
#CUTOFFVAL <- seq(140,150,20)

analpertis <- function(sampleidx, tissuedf){
    sample = SAMPLEDATA[[sampleidx]]
    print(c(tissuedf$nam[1],sample$nam[1]))
    do.call(rbind, parallel::mclapply(CUTOFFVAL,
                                      makemodel,
                                      tissue=tissuedf,
                                      sample=sample,
                                      mc.cores=12))
}

concattis <- function(tissuefilespath){
#    tissuename = getnames.unfilt(tissuefilespath)    
    tissuedf = readtissuedf(tissuefilespath)
    do.call(rbind, lapply(SAMPLE.N, analpertis, tissuedf=tissuedf))
}

runanalyses <- function(){
    tissuefiles <- list.files('/home/krishang/data/methylation/RRBS/ucscdownfiltered',pattern='.bed',full.names=T)
    tissuefiles <- tissuefiles[grepl('Osteo', tissuefiles)]
    do.call(rbind, lapply(tissuefiles, concattis))
}

#####

# concattis <- function(sampleidx){
#     sample <- SAMPLEDATA[[sampleidx]]
#     ## do.call(rbind, lapply(TISSUE.N, analpertis, sample=sample))
#     do.call(rbind, lapply(TISSUE.N, analpertis, sample=sample))
# }

# runanalyses <- function(){
#     do.call(rbind, lapply(SAMPLE.N, concattis))
# }
#####

samplefiles <- list.files('temp', pattern='bedcoord.txt.gz',full.names=T)
samplefiles <- samplefiles[grepl('RRBS', samplefiles)]
samplefiles <- samplefiles[grepl('Alt', samplefiles)]
SAMPLE.N <- sapply(samplefiles, getnames)
print(SAMPLE.N)

SAMPLEDATA <-  parallel::mclapply(samplefiles, readsampledf, mc.cores=10)
names(SAMPLEDATA) <- SAMPLE.N
SAMPLEDATA <- SAMPLEDATA[!sapply(SAMPLEDATA,is.null)]
SAMPLE.N <- names(SAMPLEDATA)
print(names(SAMPLEDATA))
print(SAMPLE.N)

bigdf <- runanalyses()

write.table(bigdf, file='only_Alt_RRBSlinearmodel_newmethyl_withvariance_bothstrands.txt',row.names=F,col.names=T,quote=F,sep='\t')

exit()

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
