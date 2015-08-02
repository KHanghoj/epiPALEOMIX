args <- commandArgs(trailingOnly = TRUE)
smooth_val <- as.integer(args[1])
IN_path = args[2]
BEDTYPE = args[3]
OUT_path = args[4]
FOURIERBOOL = args[5]
require(ggplot2)
source('~/research/projects/epiomix/epipipe/epipaleomix/downstreamanalysis/extendsmooth.R')
merge.data <- function(uniq, df){
    keep <- df$rela_pos==uniq
    deamin = sum(df[keep,3])
    tot = sum(df[keep,4])
    deamin/tot
}

merge.data1 <- function(uniq){
    keep <- a$rela_pos==uniq
    mean(a[keep,3]/a[keep,4])
}

read.data.frames <-  function(f_in, windowsize){
    df <- read.table(f_in)
    df$bedstart <- sapply(strsplit(as.character(df[,5]), "_"), "[[", 2)
    df$rela_pos <- as.integer(df[,2])-as.integer(df$bedstart)
    uniqs <- sort(unique(df$rela_pos))
    newdf <- data.frame(x=uniqs, y=unlist(parallel::mclapply(uniqs,merge.data, df=df,mc.cores=3)))
                                        #    newdf$x <- as.integer(max(newdf$x)/2) - newdf$x
    newdf$x <- as.integer(newdf$x - as.integer(max(newdf$x)/2))

    newdf$smooth <- extendsmooth(newdf$y, windowsize)
    newdf$plotname <- makename(f_in)
    newdf
}

plotting <- function(df){
    p <- ggplot(df,aes(x,smooth))+
        geom_line()+ facet_wrap(~plotname,ncol=3,scales='free_y')+
            theme(legend.position="none")+
                labs(title=choosename(BEDTYPE),
                     y='Smoothed Methylation score (MS)',
                     x='Relative position')
    p
}



files <- list.files(IN_path, pattern='MethylMap',recursive=T,full.names=T)
files <- files[grepl(BEDTYPE, files)]
print(files)
dfs <- do.call(rbind, parallel::mclapply(files, read.data.frames, windowsize=smooth_val, mc.cores=3))
pdf(OUT_path)
print(plotting(dfs))
dev.off()

if (!is.na(FOURIERBOOL)){
fourierlist <- fourier(dfs) # from the CTCF calc
pdf(paste0(unlist(strsplit(OUT_path ,'\\.'))[1], '_fourierLEFTpart.pdf'))
print(plotting.fourier(fourierlist$left, 'left'))
dev.off()

pdf(paste0(unlist(strsplit(OUT_path ,'\\.'))[1], '_fourierRIGHTpart.pdf'))
print(plotting.fourier(fourierlist$right, 'right'))
dev.off()
}
