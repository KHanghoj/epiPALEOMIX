args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]

# require(ggplot2); require(zoo); require(TSA); require(gridExtra)
require(ggplot2); require(zoo)

plotting = function(df){
 p <- ggplot(df,aes(x,smooth,col=names,
 group=names, fill=names))+
 geom_line()+ facet_wrap(~names,ncol=3)+
 theme(legend.position="none")+
 labs(title='',
 y='Merged Methylation score (MS)',
 x='Relative position')
 p
 }

# need to add x and 0 for when data not available. add a 0
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
merge.data = function(uniq, a){
 keep <- a$rela_pos==uniq
 deamin = sum(a[keep,3])
 tot = sum(a[keep,4])
 deamin/tot
 }

merge.data1 = function(uniq){
 keep <- a$rela_pos==uniq
 mean(a[keep,3]/a[keep,4])
 }

main = function(a){
 uniqs = sort(unique(a$rela_pos))
 df = data.frame(x=uniqs, y=unlist(parallel::mclapply(uniqs,merge.data, a=a,mc.cores=3)))
 remove_from_thesides = as.integer(0.05*nrow(df))
 df=df[remove_from_thesides:(nrow(df)-remove_from_thesides),]
 df$x=seq((-nrow(df)/2)+1,nrow(df)/2,1)
 df$smooth = rollmean(df$y, smooth_val, fill='extend')
 df$names = nam
 # xx=-750:-125
 # df[(df$x %in%xx),]
 df
}

smooth_val = 50
nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]

if (file.info(IN_path)$size>55){
	a=read.table(IN_path)
	a$bedstart = sapply(strsplit(as.character(a[,5]), "_"), "[[", 2)
	a$rela_pos = as.integer(a[,2])-as.integer(a$bedstart)
	df_new = main(a)
	xlong='periodogram: %1.01f peak (bases)'
	pdf(OUT_path)
	print(plotting(df_new))
	# grid.arrange(plotting(df_new),
	# 			 plot.periodo(calc.periodo(df_new, 4), 400, xlong),
	# 			 nrow=2)
	dev.off()

} else{
	pdf(OUT_path)
	plot(1, type="n", axes=F, xlab=sprintf("NODATA in %s",nam), ylab="")
	dev.off()

}
