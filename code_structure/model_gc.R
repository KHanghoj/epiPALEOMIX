require(ggplot2,quietly = TRUE)
require(reshape2, quietly = TRUE)
require(gridExtra, quietly = TRUE)
args <- commandArgs(trailingOnly = TRUE)
input_path = args[1]
pat = args[2]
outputpath = args[3]
# output_file = file.path(getwd(), sprintf('GC_model_%s.txt',pat))
listf = list.files(input_path, pattern=pat, full.names=TRUE)
# a=read.table('out_gc_content.txt')
# a=read.table('out_gc_content_labrana.txt')
write(listf,stderr())
read.dat = function(x) {
	a=read.table(x)
	a
}
b = lapply(listf, read.dat)
names(b) = 1:length(listf)
xnames = names(b)


TVs = c()
for (i in xnames){
x = data.frame(b[[i]])

grandmean = sum(x[,3])/sum(x[,4]+1)
TV = 0.5*sum((abs(x[,3]/(x[,4]+1) - grandmean) * (x[,4]+1))/ sum(x[,4]+1) )
TVs = c(TVs,TV/grandmean)  # divide by grandmean to normalize
}
# write(TVs,stderr())

maxtv_idx = which.max(TVs)
maxtv=unique(xnames)[maxtv_idx]
df = b[[as.character(maxtv)]]
maxtv = df[1,1]
df$V3 = df$V3/max(df$V3)
df$V4 = df$V4/max(df$V4)
df$ratio = df$V3/df$V4
df$ratio[is.na(df$ratio)]=0
df$gc_content = df[,2]/df[1,1]

df_for_melt = df[,c('gc_content','V3','V4')]
colnames(df_for_melt) = c('gc_content','fragment_count', 'bin_count')
melt_df = melt(df_for_melt, id=c('gc_content'))

p1 = ggplot(melt_df, aes(gc_content,value,col=variable,group=variable)) + geom_line()+
	theme(legend.position='bottom')+labs(title=maxtv, y='count',x='GC content')
zero_pad = c(-0.1,-0.05,1,1.05)
span = 0.2
winsize=maxtv
tabsize = 3

fun = function(x){
dat = unlist(df[,5])
mean(dat[x:(x+tabsize)])
}
roll = seq(0,dim(df)[1],tabsize)
tab_means = unlist(lapply(roll,fun))

names(tab_means) = (seq(0,dim(df)[1],tabsize)+(tabsize/2))
oldpts = (as.numeric(names(tab_means)))/winsize
origpts = c(oldpts,zero_pad)
rates = c(tab_means,rep(0,length(zero_pad)))
newpts = (0:winsize)/winsize
interp = loess(rates~origpts,span = span)
predvec= predict(interp,newpts)
eps = 0.0001
predvec[predvec< max(predvec,na.rm=TRUE)*eps] = 0
# predvec[predvec< max(predvec,na.rm=TRUE)*eps] = 1
# may be more correct to do like this replacing it with 1

dat = data.frame('newpts'=newpts, 'pred_ratio'=predvec,
				 'measured_ratio'=df$ratio)
m_dat = melt(dat,id=c('newpts'))
### second plot

p2 = ggplot(m_dat, aes(newpts,value, col=variable,group=variable)) + geom_line()+
	theme(legend.position='bottom')+labs(title = sprintf('Span: %f, length: %d',span,maxtv), y='ratio', x='GC content')
###  second plot


pdf(paste(outputpath, 'GCCORRECT_modelplot.pdf',sep='/'))
grid.arrange(p1,p2,nrow=1)
invisible(dev.off())

df = cbind('count'=0:winsize,'newpts'=newpts, 'pred_ratio'=predvec, 'pred_ratio_inv'=1/predvec)
df[,4][is.infinite(df[,4])] = 0
write(df[,4],stdout(), ncolumn=1)