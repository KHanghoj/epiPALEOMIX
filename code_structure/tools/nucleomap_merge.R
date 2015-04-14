require(ggplot2); require(zoo)
remove_from_thesides = 50
windowsize = 2000
smooth = 25
nam = 'Saqqaq'
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

a=read.table('CTCF_batagai_liftover_from_hg19.txt.gz')
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
df$smooth = rollmean(df$count, smooth, fill='extend')
df$names = nam
pl = plotting(df)
pl
