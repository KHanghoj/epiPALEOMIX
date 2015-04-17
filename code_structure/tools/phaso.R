args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]
CUTOFF = as.integer(args[3])
require(ggplot2); require(TSA); require(gridExtra)

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

calc.periodo = function(df, modelt) {
	mod = lm(y~poly(x,modelt),data=df,na.action=na.omit)
	df$pred = predict(mod)
	mod_pred = df$y - df$pred
	per=periodogram(mod_pred,plot=F)
	freq = 1/per$freq
	# data.frame(x=freq,y=per$spec/max(per$spec),names=df$names)
	nam = df$names[1:length(freq)]
	data.frame(x=freq,y=per$spec,names=nam)
}

periodo.short = function(df){
df_short = df[df$x<=SHORT,]
calc.periodo(df_short, 10)
}
periodo.long = function(df){
df_long = df[df$x>SHORT,]
long = calc.periodo(df_long, 2)
}
plot.periodo = function(df, size, str) {
tophitdf<-as.integer(df[df$y==max(df$y),][1])
p1=ggplot(df, aes(x,y))+ geom_line()+
	labs(x=sprintf(str,tophitdf, SHORT),
		 y='Signal power') +
	coord_cartesian(xlim=c(0,size)) +
	theme(legend.position="none")
}
nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]

if (file.info(IN_path)$size>55){
	a=read.table(IN_path)
	df = data.frame(a)
	colnames(df) = c('x','y')
	df$names = nam

	keep = !is.na(df$y)
	df=df[keep,]
	SHORT = 125
	xshort = 'periodogram: %1.01f peak (bases) measure from 0 to %1.0f'
	xlong = 'periodogram: %1.01f peak (bases) measured from %1.0f and onward'


	pdf(OUT_path)
	grid.arrange(plotting(df[df$x>=0&df$x<=CUTOFF,]), 
				 plot.periodo(periodo.short(df), 100, xshort),
			     plot.periodo(periodo.long(df), 400, xlong),
			     nrow=3)
	dev.off()
} else{
	pdf(OUT_path)
	plot(1, type="n", axes=F, xlab=sprintf("NODATA in %s",nam), ylab="")
	dev.off()
}