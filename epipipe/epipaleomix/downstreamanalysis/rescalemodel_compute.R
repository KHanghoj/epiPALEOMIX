args <- commandArgs(trailingOnly = TRUE)
## df <- fread('zcat /disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthal_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
print(c(args, 'sizetocenter: 750 or 1500', 'inbedcoord', 'intissue', 'outputmodel'))
if(length(args)==0){stop()}

library(data.table)

fmodern <- function(bin,mdf){
    mean(mdf$V5[mdf$V5bin==bin])
}
fown <- function(bin,mdf){
    mean(mdf$methylprop[mdf$V5bin==bin])
}


df <- fread(sprintf('zcat %s',args[2]),data.table=F)
## df <- fread('zcat /disk/ginzburg/data/krishang/methyl_rawdata/bedcoord_MANYBASES_SNPsextra/AltaiNeanderthalnew_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
## tissue <- fread('/disk/ginzburg/data/krishang/methylation/RRBS/ucscdownfiltered_collapsed_new/Osteobl.bed', data.table=F)
df$nampos <- sprintf('%s_%s',df[,1],df[,2]+as.integer(args[1]))

tissue <- fread(args[3], data.table=F)
tissue$nampos <- sprintf('%s_%s',tissue[,1],tissue[,2])
df <- df[df$coverage>100,]
mdfkeep <- merge(df, tissue, by='nampos')

mdfkeep$V5bin <- with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))
par <- 1:10
rescaledf <- data.frame(modern=unlist(lapply(par,fmodern, mdf=mdfkeep)), methylprop=unlist(lapply(par,fown, mdf=mdfkeep)))
rescalemodel <- lm(modern~methylprop,data=rescaledf)
save(rescalemodel, file=args[4])


## ## load('args rescale')
## mdfkeep$methylprop <- predict(rescalemodel, data.frame(methylprop=mdfkeep$methylprop))
## m.min <- rescalemodel$model[1,1]
## m.max <- rescalemodel$model[10,1]
## mdfkeep$methylprop[mdfkeep$methylprop<m.min] <- m.min
## mdfkeep$methylprop[mdfkeep$methylprop>m.max] <- m.max
## ##
