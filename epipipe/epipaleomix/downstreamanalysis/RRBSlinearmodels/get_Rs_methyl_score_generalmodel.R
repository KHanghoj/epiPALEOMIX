make.model<- function(mdf){
    lm(methylpercent ~ methylprop*cpgread*deaminatedsites*CpGsites*deaminvariance*nodeaminvariance, na.action='na.exclude',data=mdf)
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

readtissuedf <- function(f){
    df <- fread(f, data.table=FALSE)
    #df <-  read.table(f)
    colnames(df) <- c('chrom', 'pos', 'end', 'cov', 'methylpercent')
    df$nampos <- sprintf('%s_%s', df[,1], df[,2])
    df$namtis <- getnames.unfilt(f)
    df[,4:7]  ## only take the columns we need later on.
}
readtissuedf450k <- function(f){
    df <- fread(f, data.table=FALSE)
    colnames(df) <- c('chrom', 'pos', 'end', 'methylpercent')    
    df$nampos <- sprintf('%s_%s', df[,1], df[,2])
    df$namtis <- getnames.unfilt(f)
    df[,4:6]
}


readbedcoorddata <- function(fmine, cov=100){
  sample <- fread(sprintf('zcat %s',fmine),data.table=F)
  sample <- sample[sample$coverage>=cov,]
  sample$nampos = sprintf('%s_%s_%s',sample[,1],sample[,2],sample[,3])
  sample$methylprop[sample$methylprop==0] = min(sample$methylprop[sample$methylprop>0])
  sample
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


require(data.table);require(plyr);require(reshape2)
#exit()

## this is for general genomewide model
## sample <- readsampledf('/home/krishang/data/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthal_MethylMap_RRBSk1500_bedcoord.txt.gz')
## tissue <- readtissuedf('/home/krishang/data/methylation/RRBS/ucscdownfiltered_collapsed_new/Osteobl.bed')
sample <- readsampledf('/home/krishang/data/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/Saqqaq_MethylMap_methyl450k2000_bedcoord.txt.gz')
tissue <- readtissuedf450k('/home/krishang/data/methylation/450khairsample.bed')

cov <- 200
aux <- which(sample[,5] >= cov) ## it is column 5 after removing GCCONT, as nampos is no longer first column
sampleselect <- sample[aux,]
mdf <- join(sampleselect, tissue, by='nampos',type='inner')
mod <- make.model(mdf) ###  Adjusted R-squared: 0.6763; F-statistic: 2.361e+04 on 57 and 644005 DF,  p-value: < 2.2e-16
stop()
## now we have the predictors:
## fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_PROM_bedcoord.txt.gz'
fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/Saqqaq_MethylMap_PROM_bedcoord.txt.gz'
prom <- readbedcoorddata(fmine) 

prom$predicted <- predict(mod, prom)
minpredict <- min(prom$predicted[prom$predicted>0])
prom$predicted[prom$predicted<=0] <- minpredict

quant <- quantile(prom$predicted,probs=c(0.01,0.99))
prom$predicted[prom$predicted<quant[1]] <- quant[1]
prom$predicted[prom$predicted>quant[2]] <- quant[2]


## fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/AltaiNeanderthal_MethylMap_GEBO_bedcoord.txt.gz'
fmine <- '/home/krishang/data/methyl_rawdata/bedcoord_PROMGEBO_SNPsextra/Saqqaq_MethylMap_GEBO_bedcoord.txt.gz'
gebo <- readbedcoorddata(fmine)

## fixing the windowsize since the model is fit for a window of 1400 bp
promoterwindow <- 1400
gebo$size <- with(gebo, end-pos)/promoterwindow
gebo$deaminatedsites <- gebo$deaminatedsites/gebo$size
gebo$coverage <- gebo$coverage/gebo$size
gebo$CpGsites <- gebo$CpGsites/gebo$size
gebo$cpgread <- gebo$cpgread/gebo$size

gebo$predicted <- predict(mod, gebo)
minpredict <- min(gebo$predicted[gebo$predicted>0])
gebo$predicted[gebo$predicted<=0] <- minpredict

quant = quantile(gebo$predicted,probs <- c(0.01,0.99))
gebo$predicted[gebo$predicted<quant[1]] <- quant[1]
gebo$predicted[gebo$predicted>quant[2]] <- quant[2]

## merging promoter ms and gebo ms
together = read.table('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',h=T)
together = fread('/disk/ginzburg/data/krishang/bedfiles/GEBOPROMTOGETHER.INFOFILE',data.table=F)
together <- together[,7:9]
names(prom)[names(prom)=='nampos'] = 'namposPROM'
names(gebo)[names(gebo)=='nampos'] = 'namposGEBO'
prommodselect <- merge(prom, together, by.x='namposPROM', by.y='PROMREG')
bigdf <- merge(gebo, prommodselect, by.x='namposGEBO', by.y='GEBOREG', suffixes=c('GEBO','PROM'))

bigdf$modelledratios = with(bigdf, predictedGEBO/predictedPROM)
bigdf$rawratios = with(bigdf, methylpropGEBO/methylpropPROM)


expressionsdata = fread('/disk/ginzburg/data/krishang/expressiondata/ensembletorefseq.txt',data.table=F)
d = merge(expressionsdata, bigdf, by.x='V3',by.y='ID')

# this is for bone samples
bonexpressiondata = fread('../../expressiondata/bone_affy_tsv',data.table=F)
major = merge(d,bone,by.x='V1',by.y='V3')
major = major[major$V6=='high quality',]
major = major[major$V5=='present',]
major$logmodelledratios = log(major$modelledratios)
major$lograwratios = log(major$rawratios)

major <- major[major$coveragePROM>=1000,]

## this is for hair samples
d = merge(expressionsdata, bigdf, by.x='V3',by.y='ID')
hair = fread('../../expressiondata/hair_affymetrix_GSE3419.tsv',data.table=F)
hair = hair[hair[,2]=='GSM77091',]
hair = hair[hair[,11]=='high quality',]
#hair = hair[hair[,10]=='present',]

major = merge(d,hair,by.x='V1',by.y='Gene ID')
major$logmodelledratios = log(major$modelledratios)
major$lograwratios = log(major$rawratios)
affy <- range(major[,38])
minaffy = affy[1]-0.01
maxaffy = affy[2]
jumps = maxaffy-minaffy
major$binned <- cut(major[,38], breaks = seq(minaffy, maxaffy, jumps/5),labels=1:5)


## cor.test(major$V4,major$modelledratios)
## cor.test(major$V4,major$rawratios)
ggplot(major,aes(x=binned,y=lograwratios))+geom_boxplot();dev.off()
ggplot(major,aes(x=binned,y=lograwratio))+geom_boxplot();dev.off()
