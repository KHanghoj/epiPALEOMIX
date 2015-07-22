args <- commandArgs(trailingOnly = TRUE)
## df <- fread('zcat /disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthal_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
print(c(args, 'sizetocenter', 'inbedcoord', 'intissue', 'outputbedcoord'))
if(is.null(args)){stop()}
df <- fread(sprintf('zcat %s',args[2]),data.table=F)
#df <- fread('zcat /disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthalnew_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
df$nampos <- sprintf('%s_%s',df[,1],df[,2]+750)
df$nampos <- sprintf('%s_%s',df[,1],df[,2]+args[1])
##tissue <- fread('/disk/ginzburg/data/krishang/methylation/RRBS/ucscdownfiltered_collapsed_new/Osteobl.bed', data.table=F)
tissue <- fread(args[3], data.table=F)
tissue$nampos <- sprintf('%s_%s',tissue[,1],tissue[,2])
rescaledmethyl <- function(df, tissue){
    mdf <- merge(df, tissue, by='nampos')
    mdfkeep <- mdf[mdf$coverage>100,]
    mdfkeep$V5bin <- with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))
    par <- 1:10
    m <- unlist(lapply(par,fmodern, mdf=mdfkeep))
    o <- unlist(lapply(par,fown, mdf=mdfkeep)) ## this is what we need for rescaling
    rescaledf <- data.frame(modern=m, methylprop=o)
    print(rescaledf)
    rescalemodel <- lm(modern~methylprop,data=rescaledf)
    df$methylprop <- predict(rescalemodel, data.frame(methylprop=df$methylprop))
    df$methylprop[df$methylprop<m[1]] <- m[1]
    df$methylprop[df$methylprop>m[10]] <- m[10]
    df
}
outputtable <- rescaledmethyl(df,tissue)
gz1 <- gzfile("/disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthalnewrescaled_MethylMap_RRBSk1500_bedcoord.txt.gz", "w")
write.table(outputtable[,1:11], file=gz1,
            row.names=F,col.names=T,quote=F,sep='\t')
close(gz1)


## err.func <- function(X1, X2){
##     sum(sqrt((X1-X2)**2))/length(X1)
## }
## library(data.table)
## df <- fread('zcat /disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthal_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
## df <- fread('zcat /disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthalnew_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
## df$nampos <- sprintf('%s_%s',df[,1],df[,2]+750)
## ## tissue <- fread('../methylation/RRBS/onlyRRBSoverlappingPROM.bed', data.table=F)
## ## tissue <- fread('../methylation/RRBS/onlyRRBSoverlappingGEBO.bed', data.table=F)
## tissue <- fread('/disk/ginzburg/data/krishang/methylation/RRBS/ucscdownfiltered_collapsed_new/Osteobl.bed', data.table=F)
## tissue$nampos <- sprintf('%s_%s',tissue[,1],tissue[,2])

## mdf <- merge(df, tissue, by='nampos')
## mdfkeep <- mdf[mdf$coverage>100,]
## mdfkeep <- mdf[mdf$coverage>1000,]  ## this gives a more clear correlation between modern and ancient
## mdfkeep <- mdf[mdf$coverage>2000,] 

## mdfkeep$V5bin <- with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))

## rescaledmethyl <- function(df, tissue){
##     mdf <- merge(df, tissue, by='nampos')
##     mdfkeep <- mdf[mdf$coverage>100,]
##     mdfkeep$V5bin <- with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))
##     par <- 1:10
##     m <- unlist(lapply(par,fmodern, mdf=mdfkeep))
##     o <- unlist(lapply(par,fown, mdf=mdfkeep)) ## this is what we need for rescaling
##     rescaledf <- data.frame(modern=m, methylprop=o)
##     print(rescaledf)
##     rescalemodel <- lm(modern~methylprop,data=rescaledf)
##     df$methylprop <- predict(rescalemodel, data.frame(methylprop=df$methylprop))
##     df$methylprop[df$methylprop<m[1]] <- m[1]
##     df$methylprop[df$methylprop>m[10]] <- m[10]
##     df
##     ## mdfkeep$methylprop <- predict(rescalemodel, data.frame(methylprop=mdfkeep$methylprop))
##     ## mdfkeep$methylprop[mdfkeep$methylprop<m[1]] <- m[1]
##     ## mdfkeep$methylprop[mdfkeep$methylprop>m[10]] <- m[10]
##     ## mdfkeep   
## }
## outputtable <- rescaledmethyl(df,tissue)
## gz1 <- gzfile("/disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthalnewrescaled_MethylMap_RRBSk1500_bedcoord.txt.gz", "w")
## write.table(outputtable[,1:11], file=gz1,
##             row.names=F,col.names=T,quote=F,sep='\t')
## close(gz1)


## fmodern <- function(bin,mdf){
##     mean(mdf$V5[mdf$V5bin==bin])
## }
## par <- 1:10
## m <- unlist(lapply(par,fmodern, mdf=mdfkeep))

## fown <- function(bin,mdf){
##     mean(mdf$methylprop[mdf$V5bin==bin])
## }
## o <- unlist(lapply(par,fown, mdf=mdfkeep)) ## this is what we need for rescaling
## rescaledf <- data.frame(modern=m, methylprop=o)
## plot(rescaledf, xlab='Moderndata 10bins', ylab='Ms score from same 10 bins')
## dev.off()
## rescalemodel <- lm(modern~methylprop,data=rescaledf)
## ## this piece of code goes in big analyses script to rescale
## ## load('args rescale')
## mdfkeep$methylprop <- predict(rescalemodel, data.frame(methylprop=mdfkeep$methylprop))
## m.min <- rescalemodel$model[1,1]
## m.max <- rescalemodel$model[10,1]
## mdfkeep$methylprop[mdfkeep$methylprop<m.min] <- m.min
## mdfkeep$methylprop[mdfkeep$methylprop>m.max] <- m.max
## ##

## mdfkeep$methylrescaled <- predict(rescalemodel, data.frame(methylprop=mdfkeep$methylprop))
## mdfkeep$methylrescaled[mdfkeep$methylrescaled<m[1]] <- m[1]
## mdfkeep$methylrescaled[mdfkeep$methylrescaled>m[10]] <- m[10]



## summary(lm(V5~methylprop*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## dim(mdfkeep)
## summary(lm(V5~methylrescaled*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared



## summary(lm(V5~methylprop*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## [1] 0.7077669
## dim(mdfkeep)
## [1] 941800     18
## summary(lm(V5~methylrescaled*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## [1] 0.7211457

## > summary(lm(V5~methylprop*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## [1] 0.6819824
## > dim(mdfkeep)
## [1] 637049     18
## > summary(lm(V5~methylrescaled*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## [1] 0.7093277


## > summary(lm(V5~methylprop*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## [1] 0.7578274
## > dim(mdfkeep)
## [1] 25757    18
## > summary(lm(V5~methylrescaled*deaminatedsites*coverage*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep))$r.squared
## [1] 0.8014873

## mdfkeep$moderndeaminated <- round(with(mdfkeep, V4*V5),0)
## mdfkeep$methylpropbin = cut(mdfkeep$methylprop, breaks= c(-1, seq(0.0007,0.007,0.007/10)), labels=1:10)
## unlist(lapply(par,f, mdfkeep=mdfkeep))



## val <- mean(mdfkeep$V5[mdfkeep$V5bin==1])
## mdfkeep$methylprop[mdfkeep$V5<val] <- val
## val <- mean(mdfkeep$V5[mdfkeep$V5bin==10])
## mdfkeep$methylprop[mdfkeep$V5>val] <- val


## mdfkeep$methylprop[mdfkeep$methylprop>0.02] <- 0.02

## val <- mean(mdfkeep$methylprop[mdfkeep$V5bin==1])
## mdfkeep$methylprop[mdfkeep$methylprop<val] <- val
## val <- mean(mdfkeep$methylprop[mdfkeep$V5bin==10])
## mdfkeep$methylprop[mdfkeep$methylprop>val] <- val


## with(mdfkeep[1:10000,], plot(methylprop, sort(V5$bin),col='green'))
## dev.off()







## glm.out.raw = glm(cbind(moderndeaminated, V4-moderndeaminated)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,family=binomial(logit),data=mdfkeep[1:100000,])
## err.func(glm.out.raw$fitted.values, with(mdfkeep[1:100000,],moderndeaminated/V4))



## glm.out.raw.lm = lm((moderndeaminated/V4)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep[1:100000,])
## err.func(glm.out.raw.lm$fitted.values, with(mdfkeep[1:100000,],moderndeaminated/V4))
## mdfkeep$V5bin= with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))

## anova(glm.out.raw,test="Chisq")



## glm.out.raw = glm(cbind(moderndeaminated, V4-moderndeaminated)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,family=binomial(logit),data=mdfkeep[1:10000,])
## glm.out.raw.lm = lm((moderndeaminated/V4)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep[1:10000,])
## glm.out.raw.simple = lm((moderndeaminated/V4)~methylprop,data=mdfkeep[1:10000,])

## d = sort(glm.out.raw$fitted.values)
## dsimple = sort(glm.out.raw.simple$fitted.values)
## dlm = sort(glm.out.raw.lm$fitted.values)

## with(mdfkeep[1:10000,], plot(d, sort(V5),col='green'))
## with(mdfkeep[1:10000,], points(dsimple, sort(V5),col='blue'))
## with(mdfkeep[1:10000,], points(dlm, sort(V5), col='red'))
## with(mdfkeep[1:10000,], points(sort(methylprop)/max(methylprop), sort(V5)))
## dev.off()

## d = (glm.out.raw$fitted.values-mean(glm.out.raw$fitted.values))/sd(glm.out.raw$fitted.values)
## l = (glm.out.raw.lm$fitted.values-mean(glm.out.raw.lm$fitted.values))/sd(glm.out.raw.lm$fitted.values)
## real = (mdfkeep[1:10000,'V5']-mean(mdfkeep[1:10000,'V5']))/sd(mdfkeep[1:10000,'V5'])
## raw = (mdfkeep$methylprop-mean(mdfkeep$methylprop))/sd(mdfkeep$methylprop)

## plot(d,real, col='green')
## points(l,real,col='red')
## points(real,real,col='black')
## points(raw[1:10000],real,col='blue')
## dev.off()
