args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
OUT_path = args[2]

require(ggplot2); require(zoo); require(TSA); require(gridExtra)


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


smooth_val = 50
nam = unlist(strsplit(IN_path,'/'))
nam = nam[length(nam)]
nam = unlist(strsplit(nam,'\\.'))[1]

a=read.table(IN_path)
a$bedstart = sapply(strsplit(as.character(a[,5]), "_"), "[[", 2)
a$rela_pos = as.integer(a[,2])-as.integer(a$bedstart)
uniqs = sort(unique(a$rela_pos))

merge.data = function(uniq){
 keep <- a$rela_pos==uniq
 deamin = sum(a[keep,3])
 tot = sum(a[keep,4])
 deamin/tot
 }

merge.data1 = function(uniq){
 keep <- a$rela_pos==uniq
 mean(a[keep,3]/a[keep,4])
 }

main = function(uniqs, func){
 df = data.frame(x=uniqs, y=unlist(parallel::mclapply(uniqs,func,mc.cores=3)))
 remove_from_thesides = as.integer(0.05*nrow(df))
 df=df[remove_from_thesides:(nrow(df)-remove_from_thesides),]
# THIS IS CENTERING X-AXIS:
# df$x=seq((-nrow(df)/2)+1,nrow(df)/2,1)
 df$smooth = rollmean(df$y, smooth_val, fill='extend')
 df$names = nam
 df
}
df_new = main(uniqs,merge.data)
xlong='periodogram: %1.01f peak (bases)'
# print(df_new)
# pl2 = main(uniqs,merge.data1)
# require(gridExtra)
pdf(OUT_path)
grid.arrange(plotting(df_new),
			 plot.periodo(calc.periodo(df_new, 4), 400, xlong),
			 nrow=2)
# grid.arrange(pl1,pl2)
dev.off()