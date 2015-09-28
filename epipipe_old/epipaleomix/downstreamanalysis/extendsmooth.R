smoothfunc <- function(idx, seq, windowsize){
sum(seq[idx:(idx+windowsize-1)])/(windowsize)
}

extendsmooth <- function(seq, windowsize=2){
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


makename <-  function(f_in){
    nam = unlist(strsplit(as.character(f_in),'/'))
    nam = nam[length(nam)]
    as.character(unlist(strsplit(nam,'_'))[1])
}

CHANGENAMES <- list('EtoI' = 'Exon to Intron',
                    'ItoE' = 'Intron to Exon',
                    'PROM' = 'Promoters',
                    'third' = 'CpG Islands',
                    'PROMHIGH' = 'High',
                    'PROMHOUSEKEEPGENES' = 'House-Keeping',
                    'PROMINTERMEDIATE' = 'Intermediate',
                    'PROMLOW'= 'Low',
                    'CpGISLmidthird' = 'Intermediate',
                    'CpGISLtopthird' = 'Low',
                    'CpGISLlowestthird' = 'High',
                    'ISL' = 'CpG Islands with Shores (2kb) and Shelves (2kb)')


choosename <- function(n){
    if (is.null(CHANGENAMES[[n]])){n}
    else {CHANGENAMES[[n]]}
}


require(TSA)

calc.periodo <- function(n, bigdf, modelt) {
                                        # mod = lm(y~poly(x,modelt),data=df,na.action=na.omit)
    df <- bigdf[[n]]
    #mod <- lm(y~poly(x,modelt),data=df,na.action=na.omit)
    mod <- lm(smooth~poly(x,modelt),data=df,na.action=na.omit)
    df$pred <- predict(mod)
    mod_pred <- df$smooth - df$pred
    per <- periodogram(mod_pred,plot=F)
    newdf <- data.frame(x=1/per$freq, y=per$spec)
    newdf$plotname <- n
    newdf
}

plotting.fourier <- function(df, side){
    sidesinfo <- list('left' = 'from -950 to -50 relative to ',
                      'right' = 'from 50 to 950 relative to')
    p <- ggplot(df,aes(x,y)) + geom_line()+ facet_wrap(~plotname,ncol=3,scales='free_y') +
        theme(legend.position="none") +
            labs(title=sprintf('Fourier Transformation: %s %s',sidesinfo[[side]],choosename(BEDTYPE)),
                 y='Signal Power)',
                 x='WaveLength (bp)')
    p
}

splitdfs <- function(df, firstval, secondval){
    sidedf <- df[(df$x < firstval & df$x > secondval), ]
    bigdfsplit <-  split(sidedf, sidedf$plotname)
    do.call(rbind, lapply(names(bigdfsplit), calc.periodo, bigdf=bigdfsplit, modelt=4))
    
}

fourier <- function(dfs){
    list('left'=splitdfs(dfs, -50, -950),
         'right'=splitdfs(dfs, 950, 50))
}
