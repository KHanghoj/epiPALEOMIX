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
 geom_line()+ facet_wrap(~names,ncol=3)+
 theme(legend.position="none")+
 labs(title='',
 y='Merged Methylation score (MS)',
 x='Relative position')
 p
 }

# need to add x and 0 for when data not available. add a 0

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
 df = df[df$x>=0&df$x<=windowsize,]
 df$x = df$x-(windowsize/2)
 df=df[remove_from_thesides:(nrow(df)-remove_from_thesides),]
 df$smooth = rollmean(df$y, smooth_val, fill='extend')
 df$names = nam
 plotting(df)
}
pl1 = main(uniqs,merge.data)
# pl2 = main(uniqs,merge.data1)
# require(gridExtra)
pdf(OUT_path)
pl1
# grid.arrange(pl1,pl2)
dev.off()
