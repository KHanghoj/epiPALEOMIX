args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]

require(ggplot2); require(TSA); require(zoo); require(gridExtra)
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

calc.periodo = function(df, modelt) {
	# mod = lm(y~poly(x,modelt),data=df,na.action=na.omit)
	mod = lm(smooth~poly(x,modelt),data=df,na.action=na.omit)
	df$pred = predict(mod)
	mod_pred = df$smooth - df$pred
	per=periodogram(mod_pred,plot=F)
	freq = 1/per$freq
	# data.frame(x=freq,y=per$spec/max(per$spec),names=df$names)
	nam = df$names[1:length(freq)]
	data.frame(x=freq,y=per$spec,names=nam)
}
plot.periodo = function(df, size, str) {
tophitdf<-as.integer(df[df$y==max(df$y),][1])
p1=ggplot(df, aes(x,y))+ geom_line()+
	labs(x=sprintf(str,tophitdf),
		 y='Signal power') +
	coord_cartesian(xlim=c(0,size)) +
	theme(legend.position="none")
}
smooth_val = 50


nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]

if (file.info(IN_path)$size>55){
a=read.table(IN_path)
a$bedstart = sapply(strsplit(as.character(a[,6]), "_"), "[[", 2)
a$rela_pos = as.integer(a[,2])-as.integer(a$bedstart)
# if cutoff
b = a[a[,5]>=1,]
df = data.frame(table(b$rela_pos))
colnames(df) = c('x','count')
df$x = as.numeric(as.character(df$x))


remove_from_thesides = as.integer(0.05*nrow(df))
df=df[remove_from_thesides:(nrow(df)-remove_from_thesides),]
# THIS IS CENTERING X-AXIS:
# df$x=seq((-nrow(df)/2)+1,nrow(df)/2,1)

df$smooth = rollmean(df$count, smooth_val, fill='extend')
df$names = nam

xlong='periodogram: %1.01f peak (bases) measured across the entrie sample'

pdf(OUT_path)
grid.arrange(plotting(df),
			 plot.periodo(calc.periodo(df, 4), 400, xlong),
			 nrow=2)
# grid.arrange(pl1,pl2)
dev.off()
} else{
	pdf(OUT_path)
	plot(1, type="n", axes=F, xlab=sprintf("NODATA in %s",nam), ylab="")
	dev.off()
}