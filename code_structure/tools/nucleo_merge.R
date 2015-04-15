args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]

require(ggplot2); require(zoo)
remove_from_thesides = 75
windowsize = 2000
smooth_val = 50


nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]


plotting = function(df){
 p <- ggplot(df,aes(x,smooth,col=names,
 group=names, fill=names))+
 geom_area()+ facet_wrap(~names,ncol=3)+
 theme(legend.position="none")+
 labs(title='',
 y='Predicted Nucleosome Count',
 x='Relative position')
 p
}

# need to add x and 0 for when data not available. add a 0

#a=read.table('Saqqaq_CTCFbed1_NucleoMap.txt.gz')
a=read.table(IN_path)
a$bedstart = sapply(strsplit(as.character(a[,6]), "_"), "[[", 2)
a$rela_pos = as.integer(a[,2])-as.integer(a$bedstart)
# if cutoff
b = a[a[,5]>1,]
df = data.frame(table(b$rela_pos))
colnames(df) = c('x','count')
df$x = as.numeric(as.character(df$x))

df = df[df$x>=0&df$x<=windowsize,]
df$x = df$x-(windowsize/2)
df=df[remove_from_thesides:(nrow(df)-remove_from_thesides),]
df$smooth = rollmean(df$count, smooth_val, fill='extend')
df$names = nam
pl = plotting(df)
pdf(OUT_path)
pl
# grid.arrange(pl1,pl2)
dev.off()
