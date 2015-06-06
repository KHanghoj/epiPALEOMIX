args <- commandArgs(trailingOnly = TRUE)
smooth_val <- as.integer(args[1])
IN_path = args[2]
BEDTYPE = args[3]
OUT_path = args[4]

#  Rscript nucleo_merge.R 50 ~/Desktop/results/allresults/RISESAMPLEOUTPUT/ CTCF ~/Desktop/hmm.pdf

smoothfunc = function(idx, seq, windowsize){
sum(seq[idx:(idx+windowsize-1)])/(windowsize)
}

extendsmooth = function(seq, windowsize=2){
    idx = 1:length(seq)
    smooth = sapply(idx, smoothfunc, seq=seq, windowsize=windowsize)
    keep = !is.na(smooth)
    smooth = smooth[keep]
    additions = sum(!keep)
    if (additions %% 2 == 0){
        smooth = c(rep(smooth[1],additions/2),smooth, rep(smooth[length(smooth)],additions/2))
    } else {
        smooth = c(smooth, smooth[length(smooth)])
        additions = additions-1
        smooth = c(rep(smooth[1],additions/2),smooth, rep(smooth[length(smooth)],additions/2))
    }
    smooth
}

require(ggplot2) #; require(TSA); require(gridExtra)

makename <-  function(f){
    nam = unlist(strsplit(as.character(f),'/'))
    nam = nam[length(nam)]
    as.character(unlist(strsplit(nam,'_'))[1])
}

read.data.frames <-  function(file, windowsize){
    df <- read.table(file)
    df$bedstart <- sapply(strsplit(as.character(df[,6]), "_"), "[[", 2)
    df$rela_pos <- as.integer(df[,2])-as.integer(df$bedstart)
    df <- df[df[,5]>=1,]
    df <- data.frame(table(df$rela_pos))
    colnames(df) <- c('x','count')
    df$x <- as.numeric(as.character(df$x))
    df$smooth <- extendsmooth(df$count, windowsize)
    df$plotname <- makename(file)
    as.data.frame((df))
}

plotting <- function(df){
    p <- ggplot(df,aes(x,smooth))+
        geom_line()+ facet_wrap(~plotname,ncol=3,scales='free_y')+
            theme(legend.position="none")+
                labs(title=BEDTYPE,
                     y='Predicted Nucleosome Count',
                     x='Relative position')
    p
}


files = list.files(IN_path, pattern='NucleoMap',recursive=T,full.names=T)
files = files[grepl(BEDTYPE, files)]
dfs <-  do.call(rbind, lapply(files, read.data.frames, windowsize=smooth_val))
pdf(OUT_path)
plotting(dfs)
dev.off()
