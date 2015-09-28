make.model<- function(mdf){
    lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf)
}

make.big.model <- function(mdf){
    lm(m ~ methylprop*cpgreadPROM*deaminatedsitesPROM*CpGsitesPROM*deaminvariancePROM*nodeaminvariancePROM *
       cpgreadGEBO*deaminatedsitesGEBO*CpGsitesGEBO*deaminvarianceGEBO*nodeaminvarianceGEBO, na.action='na.exclude',data=mdf)
}
old.make.big.model <- function(mdf){
    lm(m ~ methylpropPROM*cpgreadPROM*deaminatedsitesPROM*CpGsitesPROM*deaminvariancePROM*nodeaminvariancePROM*methylpropGEBO*cpgreadGEBO*deaminatedsitesGEBO*CpGsitesGEBO*deaminvarianceGEBO*nodeaminvarianceGEBO, na.action='na.exclude',data=mdf)
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

f.calc <- function(idx, dfs){
       df = dfs[[idx]]
       methylreads = sum(df[,4]*df[,5])
       totalreads = sum(df[,4])
       data.frame('m'=methylreads/totalreads, 'sites'=nrow(df),'nampos'=idx)
}

readdfRRBS <- function(f){
       a=fread(f, data.table=F)
       a$nampos = with(a, sprintf('%s_%s_%s',V6,V7,V8))
       split(a,a$nampos)
}

makedf <- function(fRRBS, fmine, cov=100,promsize=1400){
  dfs = readdfRRBS(fRRBS)
  df = do.call(rbind, parallel::mclapply(names(dfs),f.calc, dfs=dfs,mc.cores=5))
  sample = fread(sprintf('zcat %s',fmine),data.table=F)
  sample = sample[sample$coverage>=cov,]
                                        # note:: takes the size of the window into account is it varies
  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
  sample$size = with(sample, end-pos)/promsize
  sample$deaminatedsites = sample$deaminatedsites/sample$size
  sample$coverage = sample$coverage/sample$size
  sample$CpGsites = sample$CpGsites/sample$size
  sample$cpgread = sample$cpgread/sample$size
  mdf = merge(sample,df,by='nampos')
  mdf
  ## mod = lm(m ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, data=mdf)
  ## print(c(getnames(fRRBS),nrow(mdf)))
  ## print(summary(mod)$r.squared)  # R-squared on .69  # R-squared on .72 housekeeping  # 
  ## mod
}


readbedcoorddata <- function(fmine, cov=100){
  sample <- fread(sprintf('zcat %s',fmine),data.table=F)
  sample <- sample[sample$coverage>=cov,]
  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
  sample
}

require(data.table);require(plyr);require(reshape2);require(ggplot2)
#exit()


## make two good predictors
fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
fRRBS <- '/disk/ginzburg/data/krishang/bedfiles/PROMRRBSOVERLAP.bed'
promdf <- makedf(fRRBS,fmine,cov=100)  ## [1] PROMRRBSOVERLAP.bed # [1] 0.7103264
## fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
## fRRBS <- '/disk/ginzburg/data/krishang/bedfiles/GEBORRBSOVERLAP.bed' 
## gebodf <- makedf(fRRBS,fmine,cov=100)
fmine='/home/krishang/data/methyl_rawdata/bedcoord_GEBOSHORT_SNPsextra/AltaiNeanderthal_MethylMap_GEBOSHORT_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBOSHORTRRBSOVERLAP.bed' 
gebodf = makedf(fRRBS,fmine,cov=100)  

together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER_GEBOSHORT.INFOFILE',h=T)
together <- together[,7:9]
names(promdf)[names(promdf)=='nampos'] = 'namposPROM'
names(gebodf)[names(gebodf)=='nampos'] = 'namposGEBO'
prommodselect <- merge(promdf, together, by.x='namposPROM', by.y='PROMREG')
bigdf <- merge(gebodf, prommodselect, by.x='namposGEBO', by.y='GEBOREG', suffixes=c('GEBO','PROM'))
bigdf$mGEBO[bigdf$mGEBO==0] = min(bigdf$mGEBO[bigdf$mGEBO!=0])
bigdf$mPROM[bigdf$mPROM==0] = min(bigdf$mPROM[bigdf$mPROM!=0])

bigdf$m = with(bigdf, mGEBO/mPROM)
bigdf$methylprop = with(bigdf, methylpropGEBO/methylpropPROM)

model <- make.big.model(bigdf)  ## 0.86 of the 2000 genes

fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
fRRBS <- '/disk/ginzburg/data/krishang/bedfiles/PROMRRBSOVERLAP.bed'
promdf <- makedf(fRRBS,fmine,cov=100)  ## [1] PROMRRBSOVERLAP.bed # [1] 0.7103264
fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
fRRBS <- '/disk/ginzburg/data/krishang/bedfiles/GEBORRBSOVERLAP.bed'
gebodf <- makedf(fRRBS,fmine,cov=100)
## fmine='/home/krishang/data/methyl_rawdata/bedcoord_GEBOSHORT_SNPsextra/AltaiNeanderthal_MethylMap_GEBOSHORT_bedcoord.txt.gz'
## fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBOSHORTRRBSOVERLAP.bed'
## gebodf = makedf(fRRBS,fmine,cov=100)

together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',h=T)
together <- together[,7:9]
names(promdf)[names(promdf)=='nampos'] = 'namposPROM'
names(gebodf)[names(gebodf)=='nampos'] = 'namposGEBO'
prommodselect <- merge(promdf, together, by.x='namposPROM', by.y='PROMREG')

bigdf <- merge(gebodf, prommodselect, by.x='namposGEBO', by.y='GEBOREG', suffixes=c('GEBO','PROM'))
bigdf$mGEBO[bigdf$mGEBO==0] = min(bigdf$mGEBO[bigdf$mGEBO!=0])
bigdf$mPROM[bigdf$mPROM==0] = min(bigdf$mPROM[bigdf$mPROM!=0])
bigdf$m = with(bigdf, mGEBO/mPROM)
bigdf$methylprop = with(bigdf, methylpropGEBO/methylpropPROM)

bigdf$predictedRATIOS <- predict(model, bigdf)
exit()

expressionsinfofile = read.table('/disk/ginzburg/data/krishang/expressiondata/ensembletorefseq.txt')
d = merge(expressionsdata, bigdf[,c('ID','predictedratio','rawratios','coveragePROM','methylpropPROM', 'predictedPROM')], by.x='V3',by.y='ID')
bonexpressiondata = fread('../../expressiondata/bone_affy_tsv',data.table=F)
bone=bonexpressiondata[bonexpressiondata$V1 =='Ost1-1',]
major = merge(d,bone,by.x='V1',by.y='V3')
major = major[major$V6=='high quality',]
#major = major[major$V5=='present',]
major$logmodelledratios = log(major$modelledratios)
major$lograwratios = log(major$rawratios)

major <- major[major$coveragePROM>=1000,]

## cor.test(major$V4,major$modelledratios)
## cor.test(major$V4,major$rawratios)

affy <- range(major$V4)
minaffy = affy[1]
maxaffy = affy[2]
jumps = maxaffy-minaffy
test <- cut(major$V4, breaks = c(-1,seq(minaffy, maxaffy, jumps/9)),labels=1:10)

major$binned <- cut(major$V4, breaks = c(-1,seq(minaffy, maxaffy, jumps/9)),labels=1:10)


major$binned <- cut(major$V4, breaks = c(-1,seq(minaffy, maxaffy, jumps/9)),labels=1:10)
major$lograwratio <- log(major$rawratio)
major$logmodelledratios <- log(major$modelledratios)


ggplot(major,aes(x=binned,y=logmodelledratios))+geom_boxplot();dev.off()
ggplot(major,aes(x=binned,y=lograwratio))+geom_boxplot();dev.off()


