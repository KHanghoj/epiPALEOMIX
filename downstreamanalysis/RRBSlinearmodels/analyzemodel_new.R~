make.model<- function(mdf){
    lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf)
}

make.big.model<- function(mdf){
    lm(methylpercent ~ methylpropPROM*cpgreadPROM*deaminatedsitesPROM*CpGsitesPROM*deaminvariancePROM*nodeaminvariancePROM*methylpropGEBO*cpgreadGEBO*deaminatedsitesGEBO*CpGsitesGEBO*deaminvarianceGEBO*nodeaminvarianceGEBO, na.action='na.exclude',data=mdf)
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
    print(nrow(df)) ## removing zeros and top values
    if(nrow(df)>5000){
        df$nampos <- sprintf('%s_%s', df[,1], df[,2]+lengthtocenter)
        df$namsam <- nam
        df$bedregion <- bedregion
        df
    }
}

# f.calc <- function(idx, dfs){
#        df = dfs[[idx]]
#        data.frame('m'=mean(df[,5]), 'v'=var(df[,5]), 'sites'=nrow(df),'nampos'=idx)
# }

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

makemodel <- function(fRRBS, fmine, cov=500){
	  dfs = readdfRRBS(fRRBS)
	  df = do.call(rbind, parallel::mclapply(names(dfs),f.calc, dfs=dfs,mc.cores=5))
	  sample = fread(sprintf('zcat %s',fmine),data.table=F)
	  sample = sample[sample$coverage>=cov,]
	  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
	  mdf = merge(sample,df,by='nampos')
	  mod = lm(m ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, data=mdf)
	  print(c(getnames(fRRBS),nrow(mdf)))
	  print(summary(mod)$r.squared)  # R-squared on .69  # R-squared on .72 housekeeping  # 
	  data.frame('fitted'=mod$fitted.values,
                     'region' = mdf$nampos,
                     'residuals'=mod$residuals,
                     'ms'=mdf$methylprop,
                     'RRBS'=mdf$m)
}

makemodel.new <- function(fRRBS, fmine, cov=100){
	  dfs = readdfRRBS(fRRBS)
	  df = do.call(rbind, parallel::mclapply(names(dfs),f.calc, dfs=dfs,mc.cores=5))
	  sample = fread(sprintf('zcat %s',fmine),data.table=F)
	  sample = sample[sample$coverage>=cov,]
	  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
	  mdf = merge(sample,df,by='nampos')
	  mod = lm(m ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, data=mdf)
	  print(c(getnames(fRRBS),nrow(mdf)))
	  print(summary(mod)$r.squared)  # R-squared on .69  # R-squared on .72 housekeeping  # 
          mod
}

makemodel.newgebo <- function(fRRBS, fmine, cov=100,promsize=1400){
	  dfs = readdfRRBS(fRRBS)
	  df = do.call(rbind, parallel::mclapply(names(dfs),f.calc, dfs=dfs,mc.cores=5))
	  sample = fread(sprintf('zcat %s',fmine),data.table=F)
	  sample = sample[sample$coverage>=cov,]
	  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
          sample$size = with(sample, end-pos)/promsize
          sample$deaminatedsites = sample$deaminatedsites/sample$size
          sample$coverage = sample$coverage/sample$size
          sample$CpGsites = sample$CpGsites/sample$size
          sample$cpgread = sample$cpgread/sample$size
	  mdf = merge(sample,df,by='nampos')
          # note:: takes the size of the window into account is it varies
	  mod = lm(m ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, data=mdf)
	  print(c(getnames(fRRBS),nrow(mdf)))
	  print(summary(mod)$r.squared)  # R-squared on .69  # R-squared on .72 housekeeping  # 
          mod
}


readbedcoorddata <- function(fmine, cov=100){
	  sample <- fread(sprintf('zcat %s',fmine),data.table=F)
	  sample <- sample[sample$coverage>=cov,]
	  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
	  sample
}

require(data.table);require(plyr)
exit()  ## stops it here

## this is for general genomewide model
sample <- readsampledf('/home/krishang/data/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthal_MethylMap_RRBSk1500_bedcoord.txt.gz')
tissue <- readtissuedf('/home/krishang/data/methylation/RRBS/ucscdownfiltered_collapsed_new/Osteobl.bed')
cov <- 1000
aux <- which(sample[,5] >= cov) ## it is column 5 after removing GCCONT, as nampos is no longer first column
sampleselect <- sample[aux,]
mdf <- join(sampleselect, tissue, by='nampos',type='inner')
mod <- make.model(mdf) ###  Adjusted R-squared: 0.6763; F-statistic: 2.361e+04 on 57 and 644005 DF,  p-value: < 2.2e-16

#### using above genome-wide model to predict ms's values in promoter and gebo regions:
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
prom = readbedcoorddata(fmine)
prom$predictedms = predict(mod, prom)
fmine='/home/krishang/data/methyl_rawdata/bedcoord_GEBOSHORT_SNPsextra/AltaiNeanderthal_MethylMap_GEBOSHORT_bedcoord.txt.gz'
gebo = readbedcoorddata(fmine)
gebo$predictedms = predict(mod, gebo)
together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER_GEBOSHORT.INFOFILE',h=T)
together <- together[,7:9]
promselect = merge(prom, together, by.x='nampos', by.y='PROMREG')
promgebotogether = merge(gebo, promselect, by.x='nampos', by.y='GEBOREG', suffixes=c('GEBO','PROM'))

require(reshape2)
df = promgebotogether[,c(12,24)]
dfm = melt(df)
ggplot(dfm,aes(value, col=variable,fill=variable))+geom_bar(binwidth=0.01,position='dodge');dev.off()

#ratio <- promgebotogether$predictedmsGEBO/promgebotogether$predictedmsPROM
#ggplot(data.frame(val=ratio[ratio>-1&ratio<250]),aes(val))+geom_bar(binwidth=1);dev.off()

## test RRBS ratios :
together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',h=T)
#together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER_GEBOSHORT.INFOFILE',h=T)
together <- together[,7:9]
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/PROMRRBSOVERLAP.bed'
dfs = readdfRRBS(fRRBS)
RRBSPROM = do.call(rbind, parallel::mclapply(names(dfs),f.calc, dfs=dfs,mc.cores=5))
#fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBOSHORTRRBSOVERLAP.bed'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBORRBSOVERLAP.bed'
dfs = readdfRRBS(fRRBS)
RRBSGEBO = do.call(rbind, parallel::mclapply(names(dfs),f.calc, dfs=dfs,mc.cores=5))
RRBSmergedp =  merge(RRBSPROM, together, by.x='nampos', by.y='PROMREG')
RRBSmerged = merge(RRBSGEBO, RRBSmergedp, by.x='nampos', by.y='GEBOREG', suffixes=c('GEBO','PROM'))

df <- RRBSmerged[,c(2,5)]
dfm <- melt(df)
ggplot(dfm,aes(value, col=variable,fill=variable))+geom_bar(binwidth=0.01,position='dodge');dev.off()
##

## starts here:  bedcoord_PROMGEBO_SNPsextra
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/PROMRRBSOVERLAP.bed'
prommod = makemodel(fRRBS,fmine,cov=100)  ## [1] PROMRRBSOVERLAP.bed # [1] 0.7103264
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBORRBSOVERLAP.bed' 
gebomod = makemodel(fRRBS,fmine,cov=100)  ## [1] GEBORRBSOVERLAP.bed  [1] 0.2457081 

together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',h=T)
together <- together[,7:9]
prommodselect = merge(prommod, together, by.x='region', by.y='PROMREG')
promgebotogether = merge(gebomod, prommodselect, by.x='region', by.y='GEBOREG', suffixes=c('GEBO','PROM'))

##### this is for short GEBO's and normal PROM
## starts here:  bedcoord_PROMGEBO_SNPsextra
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/PROMRRBSOVERLAP.bed'
prommod = makemodel(fRRBS,fmine,cov=100)  ## [1] PROMRRBSOVERLAP.bed # [1] 0.7103264
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBORRBSOVERLAP.bed' 

fmine='/home/krishang/data/methyl_rawdata/bedcoord_GEBOSHORT_SNPsextra/AltaiNeanderthal_MethylMap_GEBOSHORT_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBOSHORTRRBSOVERLAP.bed' 
gebomod = makemodel(fRRBS,fmine,cov=100)  ## [1] GEBORRBSOVERLAP.bed  [1] 0.2457081 
## ALL GEBO: [1] GEBORRBSOVERLAP.bed  [1] 0.2457081   1400 bp GEBO: [1] "GEBOSHORTRRBSOVERLAP.bed" [1] 0.6138525

together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER_GEBOSHORT.INFOFILE',h=T)
together <- together[,7:9]
prommodselect = merge(prommod, together, by.x='region', by.y='PROMREG')
promgebotogether = merge(gebomod, prommodselect, by.x='region', by.y='GEBOREG', suffixes=c('GEBO','PROM'))
#promgebotogether$ratio = with(promgebotogether, fittedGEBO/fittedPROM)
promgebotogether$ratio = with(promgebotogether, residualsGEBO/residualsPROM)
qunats = quantile(promgebotogether$ratio,probs=c(0.02,0.98))
write.table(promgebotogether[promgebotogether$ratio<=qunats[1]|promgebotogether$ratio>=qunats[2],]$ID,
	      file='IDformaltailgebovsprom.txt',row.names=F,col.names=F,quote=F,sep='\t')

## make two good predictors

fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/PROMRRBSOVERLAP.bed'
prommod.lm = makemodel.new(fRRBS,fmine,cov=100)  ## [1] PROMRRBSOVERLAP.bed # [1] 0.7103264
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBORRBSOVERLAP.bed' 
gebomod.lm = makemodel.newgebo(fRRBS,fmine,cov=100)  


## now we have the predictors:
#epredict big data:
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
prom = readbedcoorddata(fmine) 
prom$predicted = predict(prommod.lm, prom)
quant = quantile(prom$predicted,probs=c(0.01,0.99))
prom$predicted[prom$predicted<quant[1]] = quant[1]
prom$predicted[prom$predicted>quant[2]] = quant[2]
fmine='/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
gebo = readbedcoorddata(fmine)

gebo$size = with(gebo, end-pos)/1400
gebo$deaminatedsites = gebo$deaminatedsites/gebo$size
gebo$coverage = gebo$coverage/gebo$size
gebo$CpGsites = gebo$CpGsites/gebo$size
gebo$cpgread = gebo$cpgread/gebo$size

gebo$predicted = predict(gebomod.lm, gebo)
quant = quantile(gebo$predicted,probs=c(0.01,0.99))
gebo$predicted[gebo$predicted<quant[1]] = quant[1]
gebo$predicted[gebo$predicted>quant[2]] = quant[2]


## ALL GEBO: [1] GEBORRBSOVERLAP.bed  [1] 0.2457081   1400 bp GEBO: [1] "GEBOSHORTRRBSOVERLAP.bed" [1] 0.6138525
together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',h=T)
together <- together[,7:9]
prommodselect = merge(prom, together, by.x='nampos', by.y='PROMREG')
promgebotogether = merge(gebo, prommodselect, by.x='nampos', by.y='GEBOREG', suffixes=c('GEBO','PROM'))

#df = data.frame(GEBO=promgebotogether$predictedGEBO, PROM=promgebotogether$predictedPROM)
#df = data.frame(ms=promgebotogether$ratioms, modelms=promgebotogether$ratio)
#dfm <- melt(df)
#ggplot(dfm,aes(value, col=variable,fill=variable))+geom_bar(binwidth=0.01,position='dodge');dev.off()
