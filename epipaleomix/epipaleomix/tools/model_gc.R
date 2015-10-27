args <- commandArgs(trailingOnly = TRUE)
## args <- unlist(list('./', 'gc_GC', './model.txt', 'ding.pdf'))  # for testing
input_path <- args[1]
pat <- args[2]
outputpath <- args[3]
outputplot <- args[4]



listf <- list.files(input_path, pattern=pat, full.names=TRUE)
b <- lapply(listf, read.table)
names(b) <- 1:length(listf)
xnames <- names(b)


TVs <- c()
for (i in xnames){
x <- data.frame(b[[i]])
grandmean <- sum(x[,3])/sum(x[,4]+1)
TV <- 0.5*sum((abs(x[,3]/(x[,4]+1) - grandmean) * (x[,4]+1))/ sum(x[,4]+1) )
TVs <- c(TVs,TV/grandmean)  # divide by grandmean to normalize
}

df <- b[[as.character(unique(xnames)[which.max(TVs)])]]
maxtv <- df[1,1]
## df$V3 <- df$V3/sum(df$V3)  # these are not needed
## df$V4 <- df$V4/sum(df$V4)  # 
df$V3 <- df$V3/max(df$V3)
df$V4 <- df$V4/max(df$V4)
df$ratio <- df$V3/df$V4
df$ratio[is.na(df$ratio)]<-0
df$gc_content <- df[,2]/df[1,1]


pdf(outputplot)
plot(df$V2, df$V3)
lines(df$V2, df$V4)
invisible(dev.off())


zero_pad <- c(-0.1,-0.05,1,1.05)
span <- 0.2
winsize<-maxtv
tabsize <- 3

fun <- function(x, dat){
    mean(dat[x:(x+tabsize)])
}

roll <- seq(0,nrow(df),tabsize)
tabmeans <- sapply(roll,fun, dat=df[,5])
oldpts <- c((roll+(tabsize/2))/winsize, zero_pad)
rates <- c(tabmeans,rep(0,length(zero_pad)))
newpts <- (0:winsize)/winsize
interp <- loess(rates~oldpts,span = span)
predvec <- predict(interp,newpts)
eps <- 0.0001
predvec[predvec< max(predvec,na.rm=TRUE)*eps] <- 0

dat <- data.frame(
    "count"=0:winsize,
    'newpts'=newpts,
    'pred_ratio'=predvec,
    'pred_ratio_inv'=1/predvec)
upper.bound = 5 # a upper bound on the enrichment of reads

dat$pred_ratio_inv[upper.bound<dat$pred_ratio_inv] <- upper.bound  
write.table(dat, file=outputpath,row.names=F,col.names=T,quote=F,sep='\t')
