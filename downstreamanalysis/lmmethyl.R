ancient = read.table('chunksofmethylnewnew.txt')
modernb = read.table('newRRBSCHUNKS1_B.txt')
modernc = read.table('newRRBSCHUNKS1_C.txt')
modernb = modernb[modernb$V1=='chr1',]
modernc = modernc[modernc$V1=='chr1',]

modernc$V7 = sprintf('%s_%s_%s',modernc$V1,modernc$V2,modernc$V3)
modernb$V7 = sprintf('%s_%s_%s',modernb$V1,modernb$V2,modernb$V3)

modernmerged = merge(modernb,modernc,by='V7', all=T)[,c(1,7,13)]
dat = apply(modernmerged[,c(2,3)],1, function(x) mean(x,na.rm=T))
modern = data.frame('V7' = modernmerged[,1], 'V6'=dat)

forsorting = data.frame(do.call(rbind,strsplit(modernmerged$V7,'_')))
modern = modern[order(as.character(forsorting[,1]),as.numeric(as.character(forsorting[,2]))),]


ancient = ancient[ancient$V5>50,]
                                        #ancient$V7 = sprintf('%s_%s_%s',ancient$V1,ancient$V2,ancient$V3)
ancient$V7 = sprintf('chr%s_%s_%s',ancient$V1,ancient$V2,ancient$V3)

ancienttest = ancient[(ancient$V7 %in% modern$V7),]
modern = modern[modern$V7 %in% ancienttest$V7,]
summary(lm(modern$V6~ancienttest$V6))
summary(lm(modern$V6~ancienttest$V6*ancienttest$V5))

library(Hmisc)
rcorr(modern$V6, ancienttest$V6, type="pearson")




## this is for binning the data
peak = mean(ancienttest$V6)+sd(ancienttest$V6)
##bin ancient
bins = c(seq(0,peak,peak/10),1)
for (idx in 1:(length(bins)-1)){
    good = ancient$V6>=bins[idx] & ancient$V6<=bins[idx+1]
    ancienttest$V6[good] = idx-1
}
## bin modern
bins = seq(0,100,10)
moderntest = modern
for (idx in 1:(length(bins)-1)){
    good = modern$V6>=bins[idx] & modern$V6<=bins[idx+1]
    moderntest$V6[good] = idx-1
}
