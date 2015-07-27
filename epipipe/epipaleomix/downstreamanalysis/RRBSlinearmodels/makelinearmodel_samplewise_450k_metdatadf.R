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

## getnames.unfilt <- function(f){
##     nam <- unlist(strsplit(f, '/'))
##     unlist(strsplit(nam[length(nam)], '\\.'))[1]
## }

## readtissuedf <- function(f){
##     ## chrom, pos, pos, methylpercent
##     df <-  read.table(f)
##     colnames(df) <- c('chrom', 'pos', 'name', 'methylpercent')    
##     df$nampos <- sprintf('%s_%s', df[,1], df[,2])
##     df$nam <- getnames.unfilt(f)
##     print(df$nam[1])
##     df
## }

readsampledf <- function(f){
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
                          '1'=c(10,5),
                          '2'=c(20,10),
                          '3'=c(30,10),
                          '4'=c(30,20))
    
analpertis <- function(tissueidx, sample){
                                        #    tissue <- TISSUEDATA[[tissueidx]]
    tissue <- METDATA[,c(tissueidx,'nampos')]
    colnames(tissue) = c('methylpercent','nampos')
    tissue$nam <- tissueidx
    print(c(tissueidx,sample$nam[1]))
    do.call(rbind, parallel::mclapply(names(ARGUMENTLIST),
                                      makemodel,
                                      tissue=tissue,
                                      sample=sample,
                                      mc.cores=6))
}

concattis <- function(sampleidx){
    sample <- SAMPLEDATA[[sampleidx]]
    do.call(rbind, parallel::mclapply(TISSUE.N, analpertis, sample=sample, mc.cores=2))
}

runanalyses <- function(){
    do.call(rbind, lapply(SAMPLE.N, concattis))
}

samplefiles <- list.files('temp', pattern='bedcoord.txt.gz',full.names=T)
## samplefiles <- samplefiles[grepl('1500', samplefiles)]
samplefiles <- samplefiles[grepl('2000', samplefiles)]
SAMPLE.N <- sapply(samplefiles, getnames)
print(SAMPLE.N)

##GCCONT <- read.table('/home/krishang/data/methylation/methyl450data/temp/methyl450k_1500_wochr.gccontent')
GCCONT <- read.table('/home/krishang/data/methylation/methyl450data/temp/methyl450k_2000_wochr.gccontent')
colnames(GCCONT) <- c('nampos', 'GCcont', 'CpG')

SAMPLEDATA <-  parallel::mclapply(samplefiles, readsampledf, mc.cores=10)
names(SAMPLEDATA) <- SAMPLE.N
SAMPLEDATA <- SAMPLEDATA[!sapply(SAMPLEDATA,is.null)]
SAMPLE.N <- names(SAMPLEDATA)
print(names(SAMPLEDATA))
print(SAMPLE.N)

load('/home/krishang/data/methylation/methyl450data/metdataDF/metdata.Rdata')
TISSUE.N <- names(METDATA)[2:(ncol(METDATA)-3)]
## tissuefiles <- list.files('/home/krishang/data/methylation/methyl450data/downloadfiltered',pattern='.bed',full.names=T)
## #tissuefiles <- tissuefiles[grep('Gm', tissuefiles)]
## print(tissuefiles)
## TISSUE.N <- sapply(tissuefiles, getnames.unfilt)

## TISSUEDATA <-  parallel::mclapply(tissuefiles, readtissuedf, mc.cores=10)
## names(TISSUEDATA) <- TISSUE.N

mega <- runanalyses()
write.table(mega, file='pertissue450klinearmodel_2000.txt',row.names=F,col.names=T,quote=F,sep='\t')
## ## hmm whats wrong
## mega$compared <- with(mega, paste(samplename, tissuename, sep='_'))
## require(ggplot2)
## m <- subset(mega, datapoints>5000&grepl('complex|compWCpGonly',model))
## splot <- function(){
##       pdf('oneplot450k_2000.pdf')
## 	print(ggplot(m,aes(compared, rsquared, col=tissuename, shape=cov_cpg_cutoff))+geom_point() +
##       		           theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5,size=8),
##                 	   legend.position='bottom'));dev.off()
## }
## splot()
