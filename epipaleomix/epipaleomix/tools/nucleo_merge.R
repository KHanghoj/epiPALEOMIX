args <- commandArgs(trailingOnly = TRUE)
smooth_val <- as.integer(args[1])
IN_path = args[2]
OUT_path = args[3]

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

require(ggplot2)
plotting = function(df){
 p <- ggplot(df,aes(x,smooth,col=names,
 group=names, fill=names))+
 geom_line()+ facet_wrap(~names,ncol=3)+
 theme(legend.position="none")+
 labs(title='',
 y='Predicted Nucleosome Count',
 x='Relative position')
 p
}

nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]

if (file.info(IN_path)$size>55){
    a=read.table(IN_path)
    
    a$bedstart = sapply(strsplit(as.character(a[,6]), "_"), "[[", 2)
    a$rela_pos = as.integer(a[,2])-as.integer(a$bedstart)
                                        # if cutoff
    b = a[a[,5]>=0,]
    df = data.frame(table(b$rela_pos))
    colnames(df) = c('x','count')
    df$x = as.numeric(as.character(df$x))
    df$x = df$x - median(df$x)

    df$smooth = extendsmooth(df$count, smooth_val)
    df$names = nam

    pdf(OUT_path)
    print(plotting(df))
    dev.off()
} else{
    pdf(OUT_path)
    plot(1, type="n", axes=F, xlab=sprintf("NODATA in %s",nam), ylab="")
    dev.off()
}
