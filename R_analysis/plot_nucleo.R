require(ggplot2)
files = list.files(pattern='*Saq*')


for (fi in files){
pilefile = paste0('../pileups/',fi)
a = read.table(pilefile)
p = ggplot(data=a,aes(V2,V3))
map = read.table(fi)
    pdf(paste0('../plots/nucleomap',paste0(fi,'.pdf')))
    print(p + geom_line() + geom_point(data=map,aes(V2,V4),col='blue') + labs(y='depth',x='chromosome 12') + coord_cartesian(ylim = c(0, max(map[,4])+5)))
    dev.off()
print(fi)
print(dim(map))
print(mean(map$V5))
}
