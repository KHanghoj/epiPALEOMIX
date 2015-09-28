args <- commandArgs(trailingOnly = TRUE)
IN_path <- args[1]
BEDTYPE <- args[2]
OUT_path <- args[3]
require(ggplot2)
source('extendsmooth.R')

read.data.frames <-  function(file, windowsize){
    df <- read.table(file)
    c <- cbind(df,do.call(rbind,strsplit(as.character(df$V5), '_')))
    c$leftside <- as.integer(as.character(c$V2))-as.integer(as.character(c[,7]))
    c$rightside <- as.integer(as.character(c$V2))-as.integer(as.character(c[,8]))
    c <- c[!(c$V4>=5 & c$V3/c$V4==1),]
    c$V5 = as.character(c$V5)
    outerleft <- c[c$leftside <= 2000,]
    innerleft <- c[c$leftside > 2000 & c$leftside<=4000,]
    outerright <- c[c$rightside >= -2000,]
    innerright <- c[c$rightside < -2000 & c$rightside>= -4000,]
    CGI <- c[!(rownames(c) %in% c(rownames(outerleft), rownames(innerleft), rownames(innerright), rownames(outerright))), ]
    df <- data.frame(
        'y'=rbind(
            sum(outerleft$V3)/sum(outerleft$V4),
            sum(innerleft$V3)/sum(innerleft$V4),
            sum(CGI$V3)/sum(CGI$V4),
            sum(innerright$V3)/sum(innerright$V4),
            sum(outerright$V3)/sum(outerright$V4)),
        'x'=c('outerleft', 'innerleft', 'CGI', 'innerright', 'outerright'),
        'plotname' = makename(file)
        )
    df$x = with(df, factor(x, levels = c('outerleft',  'innerleft',  'CGI', 'innerright', 'outerright')))
    df
}

plotting <- function(df){
    p <- ggplot(df,aes(x,y))+geom_bar(position='dodge', stat='identity')+
        facet_wrap(~plotname,ncol=3,scales='free_y')+	
            theme(legend.position="none")+		
                labs(title=choosename(BEDTYPE),		
                     y='Methylation Score')		
    p
}


files = list.files(IN_path, pattern='MethylMap',recursive=T,full.names=T)
files = files[grepl(BEDTYPE, files)] # ISL
files = files[grepl('Saqqaq',files)]
print(files)
dfs <-  do.call(rbind, parallel::mclapply(files, read.data.frames, mc.cores=10))
pdf(OUT_path)
print(plotting(dfs))
dev.off()
