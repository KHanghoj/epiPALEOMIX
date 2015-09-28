getnames <- function(f, col=1){
    nam <- unlist(strsplit(f, '/'))
    unlist(strsplit(nam[length(nam)], '_'))[col]
}

rescale <- function(df, nam){
  load(RESCALEMODELS[grepl(getnames(nam), RESCALEMODELS)])  # loads the specifc rescalemodel
  df$methylprop <- predict(rescalemodel, data.frame(methylprop=df$methylprop))
  m.min <- rescalemodel$model[1,1]
  m.max <- rescalemodel$model[10,1]
  df$methylprop[df$methylprop<m.min] <- m.min
  df$methylprop[df$methylprop>m.max] <- m.max
  df
}

readbedcoorddata <- function(fmine, GEBOPROM, cov=100){
  sample <- fread(sprintf('zcat %s',fmine),data.table=F)
  nam <- getnames(fmine)
  sample$nam <- nam
  sample <- sample[sample$coverage>=cov,]
  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
  names(sample)[names(sample)=='nampos'] = GEBOPROM   # 'namposPROM' OR namposGEBO
  rescale(sample, nam)
}

merge.func <- function(df.prom, df.gebo){
  prommodselect <- merge(df.prom, COMBINER.TABLE, by.x='namposPROM', by.y='PROMREG')
  mdf <- merge(df.gebo, prommodselect, by.x='namposGEBO', by.y='GEBOREG', suffixes=c('GEBO','PROM'))
  mdf$RsSCORES <- with(mdf, methylpropGEBO/methylpropPROM)
  mdf$RsSCORESNORM <- scale(mdf$RsSCORES)  ## mean:0, sd=1
  mdf
}

run <- function(idx){
  print(getnames(PROMS[idx]))
  merge.func(readbedcoorddata(PROMS[idx], 'namposPROM'),   readbedcoorddata(GEBOS[idx], 'namposGEBO'))
}

require(data.table);require(plyr);require(reshape2);require(ggplot2)
##COMBINER.TABLE= fread('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',data.table=F)[,7:9]
COMBINER.TABLE= fread('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE2000',data.table=F)[,7:9]

RESCALEMODELS <- list.files('/disk/ginzburg/data/krishang/methylation/rescalemodels', pattern='.Rdata',full.names=TRUE)
PROMS <- list.files('/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra', pattern='_PROM_', full.names=TRUE)
print(PROMS)
## GEBOS <- list.files('/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra', pattern='_GEBO_', full.names=TRUE)
GEBOS <- list.files('/home/krishang/data/methyl_rawdata/bedcoord_GEBO2000_SNPsextra', pattern='_GEBO_', full.names=TRUE)
print(GEBOS)
idxes <- 1:length(PROMS)
bigdf <- do.call(rbind, lapply(idxes, run))
print(dim(bigdf))
gz1 <- gzfile("/disk/ginzburg/data/krishang/methylation/Rsscores2000.txt.gz", "w")
write.table(bigdf, file=gz1,
            row.names=F,col.names=T,quote=F,sep='\t')
close(gz1)

stop()

fPROM <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
#fPROM <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/UstIshim_MethylMap_PROM_bedcoord.txt.gz'
dfPROM <- readbedcoorddata(fPROM)
fGEBO <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
#fGEBO <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/UstIshim_MethylMap_GEBO_bedcoord.txt.gz'
fGEBO <- '~/data/methyl_rawdata/bedcoord_GEBO2000_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
dfGEBO <- readbedcoorddata(fGEBO)


fPROM <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/Saqqaq_MethylMap_PROM_bedcoord.txt.gz'
#fPROM <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/UstIshim_MethylMap_PROM_bedcoord.txt.gz'
dfPROM <- readbedcoorddata(fPROM)
fGEBO <- '/home/krishang/data/methyl_rawdata/bedcoord_GEBO2000_SNPsextra/Saqqaq_MethylMap_GEBO_bedcoord.txt.gz'
#fGEBO <- '~/data/methyl_rawdata/bedcoord_GEBOSHORT_SNPsextra/Saqqaq_MethylMap_GEBOSHORT_bedcoord.txt.gz'
#fGEBO <- '~/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/Saqqaq_MethylMap_GEBO_bedcoord.txt.gz'
dfGEBO <- readbedcoorddata(fGEBO)



#together = fread('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE2000',data.table=F)
#together = fread('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER_GEBOSHORT.INFOFILE',data.table=F)
COMBINER.TABLE= fread('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',data.table=F)
COMBINER.TABLE <- COMBINER.TABLE[,7:9]
RESCALEMODELS <- list.files('/disk/ginzburg/data/krishang/methylation/rescalemodels', pattern='.Rdata',full.names=TRUE)


                                        #together <- together[,7:9]a

    
#names(dfPROM)[names(dfPROM)=='nampos'] = 'namposPROM'
#names(dfGEBO)[names(dfGEBO)=='nampos'] = 'namposGEBO'
prommodselect <- merge(dfPROM, together, by.x='namposPROM', by.y='PROMREG')
bigdf <- merge(dfGEBO, prommodselect, by.x='namposGEBO', by.y='GEBOREG', suffixes=c('GEBO','PROM'))
#bigdf$rawratios <- with(bigdf, methylpropGEBO/methylpropPROM)
bigdf$rawratios <- with(bigdf, methylpropGEBO-methylpropPROM)

stop()
bigdf$methylpropPROM = bigdf$deaminatedsitesPROM/bigdf$coveragePROM
bigdf$methylpropGEBO = bigdf$deaminatedsitesGEBO/bigdf$coverageGEBO
bigdf$rawratios <- with(bigdf, methylpropGEBO/methylpropPROM)
bigdf <- bigdf[!(is.na(bigdf$rawratios)),]




stop()
promcat = read.table('/disk/ginzburg/data/krishang/bedfiles/PROM_autosom_wochr.promcat')
bigdfpromcat <- merge(bigdf, promcat, by.x='namposPROM',by.y='V1')
methylpropPROM
methylpropGEBO
log(ratios)
bigdfpromcat = bigdfpromcat[bigdfpromcat$rawratios<14,]
ggplot(bigdfpromcat,aes(x=rawratios)) +
  stat_density(aes(ymax = ..density..,  ymin = -..density.., fill = V2, color = V2),
               geom = "ribbon", position = "identity") +
  facet_grid(~V2)+coord_flip() +
  theme(legend.position="none", axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  labs(title='Rs score (Genebody/promoter)\n splitted by promoter categories', y='', x='Rs score');dev.off()


stop()

#bigdf <- bigdf[bigdf$coveragePROM>1000,]
## this is for bone samples
## expressionsdata = read.table('/disk/ginzburg/data/krishang/expressiondata/ensembletorefseq.txt')
expressionsdata = fread('/disk/ginzburg/data/krishang/expressiondata/ensembletorefseq.txt',data.table=F)

d = merge(expressionsdata, bigdf, by.x='V3',by.y='ID')

bonexpressiondata = fread('../../expressiondata/bone_affy_tsv',data.table=F)
bonexpressiondata = fread('../../expressiondata/brain_affy_tsv',data.table=F)
print(table(keep <- bonexpressiondata$V2 %in% fread('../../expressiondata/differentexpression1',data.table=F)[,1]))
bonexpressiondata <- bonexpressiondata[keep,]
bone = bonexpressiondata
bone=bonexpressiondata[bonexpressiondata$V1 =='Ost1-1',]
bone=bonexpressiondata[bonexpressiondata$V1 =='GSM774041_Brain_1a_A',]


major = merge(d,bone,by.x='V1',by.y='V3')
major = major[major$V6=='high quality',]
#major = major[major$V5=='present',]


## this is for hair samples
expressionsdata = fread('/disk/ginzburg/data/krishang/expressiondata/ensembletorefseq.txt',data.table=F)
d = merge(expressionsdata, bigdf, by.x='V3',by.y='ID')
hair = fread('../../expressiondata/hair_affymetrix_GSE3419.tsv',data.table=F)

hair = hair[hair[,2]=='GSM77091',]
hair = hair[hair[,11]=='high quality',]
##hair = hair[hair[,10]=='present',]

major = merge(d,hair,by.x='V1',by.y='Gene ID')


## pedersen stuff
hair = fread('../../expressiondata/pedersenetalexpressiondata/hairtranscript.txt',data.table=F)
lee.et.al = read.table('/disk/ginzburg/data/krishang/expressiondata/IPIconversion/converttoNMfromIPI.txt')


major = merge(d,hair,by.x='V3',by.y='V3')  ## for pedersen data set


##ggplot(major,aes(x=binned,y=lograwratios))+geom_boxplot();dev.off()

## fmine='/home/krishang/data/methyl_rawdata/bedcoord_GEBOSHORT_SNPsextra/AltaiNeanderthal_MethylMap_GEBOSHORT_bedcoord.txt.gz'
## fRRBS = '/disk/ginzburg/data/krishang/bedfiles/GEBOSHORTRRBSOVERLAP.bed' 
## gebomod.lm = makemodel(fRRBS,fmine,cov=100)

f <- function(bin,df, col=methylpropPROM){
  median(df[df$binned==bin,col])
}
f <- function(bin,df, col=methylpropPROM){
  mean(df[df$binned==bin,col])
}

f <- function(bin,df){
  d <- df[df$binned==bin,]
  s <- nrow(d)
  sum(d[,35] == 'present')/s
}
f <- function(bin,df){
  ### this is for bone
  d <- df[df$binned==bin,]
  s <- nrow(d)
  sum(d[,30] == 'present')/s
}


unlist(lapply(sort(unique(major$binned)), f, df=major))

# for Bgee data:
br = quantile(major[,29], probs=seq(0,1,0.2))
br[1] = br[1]-0.1
major$binned <- cut(major[,29], breaks = br,labels=1:5)
barplot(unlist(lapply(sort(unique(major$binned)), f, df=major[!is.infinite(major$rawratios),], col='rawratios')),names.arg=c(1:5));dev.off()


br = quantile(major[,34], probs=seq(0,1,0.1))
br[1] = br[1]-0.1
major$binned <- cut(major[,34], breaks = br,labels=1:10)
barplot(unlist(lapply(sort(unique(major$binned)), f, df=major[!is.infinite(major$rawratios),], col='rawratios')),names.arg=c(1:10));dev.off()


# test of present vs absent
br = quantile(major$rawratios, probs=seq(0,1,0.2))
br[1] = br[1]-0.1
major$binned <- cut(major$rawratios, breaks = br,labels=1:5)
barplot(unlist(lapply(sort(unique(major$binned)), f, df=major)),names.arg=c(1:5));dev.off()

# for pederesenetaldata:
br = quantile(major$V2.y, probs=seq(0,1,0.05))
br[1] = br[1]-0.1
major$binned <- cut(major$V2.y, breaks = br,labels=1:20)
barplot(unlist(lapply(sort(unique(major$binned)), f, df=major[!is.infinite(major$rawratios),], col='rawratios')),names.arg=c(1:20));dev.off()

br = quantile(major$V2.y, probs=seq(0,1,0.2))
br[1] = br[1]-0.1
major$binned <- cut(major$V2.y, breaks = br,labels=1:5)
barplot(unlist(lapply(sort(unique(major$binned)), f, df=major[!is.infinite(major$rawratios),], col='rawratios')),names.arg=c(1:5));dev.off()


br = quantile(major$V2.y, probs=seq(0,1,0.1))
br[1] = br[1]-0.1
major$binned <- cut(major$V2.y, breaks = br,labels=1:10)
barplot(unlist(lapply(sort(unique(major$binned)), f, df=major[!is.infinite(major$rawratios),], col='rawratios')),names.arg=c(1:10));dev.off()


ggplot(major,aes(x=binned,y=rawratios))+geom_boxplot();dev.off()

barplot(unlist(lapply(sort(unique(major$binned)), f, df=major, col='rawratios')),names.arg=c(1:5));dev.off()



