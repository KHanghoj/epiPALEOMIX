args <- commandArgs(trailingOnly = TRUE)
## args <- unlist(list(c('./gc', 'gc_GC', './model.txt', 'ding.pdf')))  # for testing
input_path <- args[1]
pat <- args[2]
outputpath <- args[3]
outputplot <- args[4]


smoothfunc <- function(idx, seq, windowsize){
sum(seq[idx:(idx+windowsize-1)])/(windowsize)
}

extendsmooth <- function(seq, windowsize=2){
    idx <- 1:length(seq)
    smooth <- sapply(idx, smoothfunc, seq=seq, windowsize=windowsize)
    keep <- !is.na(smooth)
    smooth <- smooth[keep]
    additions <- sum(!keep)
    if (additions %% 2 == 0){
        smooth <- c(rep(smooth[1],additions/2),smooth, rep(smooth[length(smooth)],additions/2))
    } else {
        smooth <- c(smooth, smooth[length(smooth)])
        additions <- additions-1
        smooth <- c(rep(smooth[1],additions/2),smooth, rep(smooth[length(smooth)],additions/2))
    }
    smooth
}


listf <- list.files(input_path, pattern=pat, full.names=TRUE)
print(listf)
b <- lapply(listf, read.table)
names(b) <- 1:length(listf)
xnames <- names(b)


TVs <- c()
for (i in xnames){
x <- data.frame(b[[i]])
grandmean <- sum(x[,3]+1)/sum(x[,4]+1)
TV <- 0.5*sum((abs(x[,3]/(x[,4]+1) - grandmean) * (x[,4]+1))/ sum(x[,4]+1) )
TVs <- c(TVs,TV/grandmean)  # divide by grandmean to normalize
}

df <- b[[as.character(unique(xnames)[which.max(TVs)])]]
maxtv <- df[1,1]
colnames(df) <- c('gc', 'gccont', 'reads', 'ref')
df$readsdensity <- with(df,reads/sum(reads))+0.001
df$refdensity <- with(df,ref/sum(ref))+0.001

## this is new
df$readsdensity <- extendsmooth(df$readsdensity,3)
df$refdensity <- extendsmooth(df$refdensity,3)


df$ratio = with(df, refdensity/readsdensity)
df$gc_content <- df[,2]/df[1,1]

yranges <- c(0,with(df,max(readsdensity,refdensity)) + 0.01)
xranges <- c(0,1)

pdf(outputplot)
par(mfrow=c(1,2))
plot(1,xlim=xranges, ylim=yranges, ylab='Density', xlab='% GC content', main=sprintf('Max TV with readlength: %s bp', maxtv))
with(df, lines(gc_content,readsdensity, col='red'))
with(df, lines(gc_content,refdensity, col='blue'))
legend('topright', c('ReadDensity', 'RefDensity') ,
          lty=1, col=c('red', 'blue'), bty='n', cex=.75)

upper.bound = 5 # a upper bound on the enrichment of reads
lower.bound = 1/upper.bound
df$ratio[df$ratio > upper.bound] = upper.bound
df$ratio[df$ratio < lower.bound] = lower.bound

with(df, plot(gc_content, ratio, col='black', main='Model', ylab='correction value', xlab='% GC content',ylim=c(0,upper.bound+0.2)))
abline(1,0)
invisible(dev.off())

dat <- data.frame(
    "count"=0:maxtv,
    'GC_content'=df$gc_content,
    'pred_ratio'=1/df$ratio,
    'pred_ratio_inv'=df$ratio)

write.table(dat, file=outputpath,row.names=F,col.names=T,quote=F,sep='\t')
