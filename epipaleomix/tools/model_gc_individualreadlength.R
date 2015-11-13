args <- commandArgs(trailingOnly = TRUE)
## args <- unlist(list(c('./OUT_mk_yaml2/TEMPORARYFILES_mk_yaml2/BAMName1/', '_GC', './model.txt', 'ding.pdf')))  # for testing
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

func <- function(idx, data){
df <- data[[idx]]
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


par(mfrow=c(1,2))
plot(1,xlim=xranges, ylim=yranges, ylab='Density', xlab='% GC content', main=sprintf('Max TV with readlength: %s bp', maxtv))
with(df, lines(gc_content,readsdensity, col='red'))
with(df, lines(gc_content,refdensity, col='blue'))
legend('topright', c('ReadDensity', 'RefDensity') ,
          lty=1, col=c('red', 'blue'), bty='n', cex=.75)

upper.bound = 4 # a upper bound on the enrichment of reads
lower.bound = 1/upper.bound
df$ratio[df$ratio > upper.bound] = upper.bound
df$ratio[df$ratio < lower.bound] = lower.bound

with(df, plot(gc_content, ratio, col='black', main='Model', ylab='correction value', xlab='% GC content',ylim=c(0,upper.bound+0.2)))
abline(1,0)

dat <- data.frame(
    'count'=maxtv,
    'GC_content'=df$gc_content,
    'ratio'=df$ratio)
dat
}

listf <- list.files(input_path, pattern=pat, full.names=TRUE)
DATA <- lapply(listf, read.table)
names(DATA) <- 1:length(listf)
print(listf)
pdf(outputplot)
dat <- do.call(rbind,lapply(names(DATA), func, data=DATA))
invisible(dev.off())


write.table(dat, file=outputpath,row.names=F,col.names=T,quote=F,sep='\t')
