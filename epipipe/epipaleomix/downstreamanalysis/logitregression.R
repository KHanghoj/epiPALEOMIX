err.func <- function(X1, X2){
    sum(sqrt((X1-X2)**2))/length(X1)
}
library(data.table)
df <- fread('zcat bedcoord_MANYBASES_SNPsextra/AltaiNeanderthal_MethylMap_RRBSk1500_bedcoord.txt.gz',data.table=F)
df$nampos <- sprintf('%s_%s',df[,1],df[,2]+750)
tissue <- fread('../methylation/RRBS/onlyRRBSoverlappingPROM.bed', data.table=F)
#tissue <- fread('../methylation/RRBS/onlyRRBSoverlappingGEBO.bed', data.table=F)

tissue$nampos <- sprintf('%s_%s',tissue[,1],tissue[,2])

mdf <- merge(df, tissue, by='nampos')
mdfkeep = mdf[mdf$coverage>100,]
mdfkeep$moderndeaminated = round(with(mdfkeep, V4*V5),0)
mdfkeep$V5bin= with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))

val <- mean(mdfkeep$methylprop[mdfkeep$V5bin==1])
mdfkeep$methylprop[mdfkeep$methylprop<val] <- val
val <- mean(mdfkeep$methylprop[mdfkeep$V5bin==10])
mdfkeep$methylprop[mdfkeep$methylprop>val] <- val


with(mdfkeep[1:10000,], plot(methylprop, sort(V5$bin),col='green'))
dev.off()







glm.out.raw = glm(cbind(moderndeaminated, V4-moderndeaminated)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,family=binomial(logit),data=mdfkeep[1:100000,])
err.func(glm.out.raw$fitted.values, with(mdfkeep[1:100000,],moderndeaminated/V4))



glm.out.raw.lm = lm((moderndeaminated/V4)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep[1:100000,])
err.func(glm.out.raw.lm$fitted.values, with(mdfkeep[1:100000,],moderndeaminated/V4))
mdfkeep$V5bin= with(mdfkeep, cut(V5, breaks=c(-1,seq(0.1,1,0.1)),labels=1:10))

anova(glm.out.raw,test="Chisq")



glm.out.raw = glm(cbind(moderndeaminated, V4-moderndeaminated)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,family=binomial(logit),data=mdfkeep[1:10000,])
glm.out.raw.lm = lm((moderndeaminated/V4)~methylprop*deaminvariance*CpGsites*deaminvariance*nodeaminvariance,data=mdfkeep[1:10000,])
glm.out.raw.simple = lm((moderndeaminated/V4)~methylprop,data=mdfkeep[1:10000,])

d = sort(glm.out.raw$fitted.values)
dsimple = sort(glm.out.raw.simple$fitted.values)
dlm = sort(glm.out.raw.lm$fitted.values)

with(mdfkeep[1:10000,], plot(d, sort(V5),col='green'))
with(mdfkeep[1:10000,], points(dsimple, sort(V5),col='blue'))
with(mdfkeep[1:10000,], points(dlm, sort(V5), col='red'))
with(mdfkeep[1:10000,], points(sort(methylprop)/max(methylprop), sort(V5)))
dev.off()

d = (glm.out.raw$fitted.values-mean(glm.out.raw$fitted.values))/sd(glm.out.raw$fitted.values)
l = (glm.out.raw.lm$fitted.values-mean(glm.out.raw.lm$fitted.values))/sd(glm.out.raw.lm$fitted.values)
real = (mdfkeep[1:10000,'V5']-mean(mdfkeep[1:10000,'V5']))/sd(mdfkeep[1:10000,'V5'])
raw = (mdfkeep$methylprop-mean(mdfkeep$methylprop))/sd(mdfkeep$methylprop)

plot(d,real, col='green')
points(l,real,col='red')
points(real,real,col='black')
points(raw[1:10000],real,col='blue')
dev.off()
