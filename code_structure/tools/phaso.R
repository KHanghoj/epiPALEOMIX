args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]
CUTOFF = as.integer(args[3])
require(ggplot2); require(zoo)
nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]

plotting = function(df){
 p <- ggplot(df,aes(x,y,col=names,
 group=names, fill=names))+
 geom_line()+ facet_wrap(~names,ncol=3)+
 theme(legend.position="none")+
 labs(title='',
 y='Length Count',
 x='Relative position')
 p
}

a=read.table(IN_path)
df = data.frame(a)
colnames(df) = c('x','y')
df = df[df$x>=0&df$x<=CUTOFF,]
df$names = nam
pl = plotting(df)
pdf(OUT_path)
pl
# grid.arrange(pl1,pl2)
dev.off()
