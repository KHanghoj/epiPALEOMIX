f.mostcomplex <- function(mdf){
    summary(lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf))$r.squared
}

getnames <- function(f, col=1){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[col]
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
    bedregion <- getnames(f, col=3)
    print(nam)
    df <- fread(sprintf('zcat %s', f), data.table=FALSE,h=T)
    load(RESCALEMODELS[grepl(getnames(nam), RESCALEMODELS)])  # loads the specifc rescalemodel
    df$methylprop <- predict(rescalemodel, data.frame(methylprop=df$methylprop))
    m.min <- rescalemodel$model[1,1]
    m.max <- rescalemodel$model[10,1]
    df$methylprop[df$methylprop<m.min] <- m.min
    df$methylprop[df$methylprop>m.max] <- m.max
    print(nrow(df)) ## removing zeros and top values
    df$nampos <- sprintf('%s_%s', df[,1], df[,2]+lengthtocenter)
    df$namsam <- nam
    df$bedregion <- bedregion
    df
}

makemodel <- function(cov, tissue, sample){
    aux <- which(sample[,5] >= cov) ## it is column 5 after removing GCCONT, as nampos is no longer first column
    sampleselect <- sample[aux,]
    print(c(tissue$namtis[1],sample$namsam[1], cov))
    mdf <- join(sampleselect, tissue, by='nampos',type='inner')
    if(nrow(mdf)>100){data.frame('model'=c('Mostcomplex'),
                                 'rsquared'=f.mostcomplex(mdf),
                                 'datapoints'=nrow(mdf),
                                 'cov_cutoff'=cov,
                                 'samplename'=mdf$namsam[1],
                                 'tissuename'=mdf$namtis[1],
                                 'bedregion'=mdf$bedregion[1]
                                 )
                  }
}


analpertis <- function(tissuefilespath, sample){
    tissuedf <- readtissuedf(tissuefilespath)
    forappend <- do.call(rbind, parallel::mclapply(CUTOFFVAL,
                                                   makemodel,
                                                   tissue=tissuedf,
                                                   sample=sample,
                                                   mc.cores=10))
    write.table(forappend, file='ALLSAMPLES_COMPLEX_RRBSlinearmodel_MANYBASES_MERGEDRRBSDATA_ALLTISSUES_RESCALE_SNPREMOVE.txt',append=T,row.names=F,col.names=F,quote=F,sep='\t')
}

runanalyses<- function(samplepath){
    sample <- readsampledf(samplepath)
    tissuefiles <- list.files('/home/krishang/data/methylation/RRBS/ucscdownfiltered_collapsed',pattern='.bed',full.names=T)
    tissuefiles <- tissuefiles[grepl('.bed$',tissuefiles)]
    #tissuefiles <- tissuefiles[grepl('Osteobl', tissuefiles)]
    parallel::mclapply(tissuefiles, analpertis, sample=sample,mc.cores=6)
}

require(data.table)
require(plyr)

RESCALEMODELS <- list.files('/disk/ginzburg/data/krishang/methylation/rescalemodels', pattern='.Rdata',full.names=TRUE)
#CUTOFFVAL <- c(seq(1,2500,25))
CUTOFFVAL <- c(seq(100,2500,50))
samplefiles <- list.files('/home/krishang/data/methyl_rawdata/bedcoord_MANYBASES_SNPsextra', pattern='bedcoord.txt.gz',full.names=T)
samplefiles <- samplefiles[grepl('RRBS', samplefiles)]
samplefiles <- samplefiles[grepl('1500', samplefiles)]
#samplefiles <- samplefiles[grepl('Altai', samplefiles)] 
print(samplefiles)
mega <- do.call(rbind, lapply(samplefiles, runanalyses))

quit('no')
