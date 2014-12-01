# require(ggplot2)
# require(zoo)
files = list.files(pattern='*Saq*')
Test <- lapply(files[1], function(x) {
  a <- read.table(x)
  start_nucleo = a$V2
  i = 1
  for (idx in 1:length(start_nucleo)){
  	if (i > 1){  # i do not want to get the same hit several times. so pass loop == the length of former nuclesome phasing
  		i = i-1
  		next
  	}
  	currentstart = start_nucleo[idx]
  	while (start_nucleo[idx+i] <= currentstart+1000 & start_nucleo[idx+i] >= start_nucleo[idx+i-1]+147 & i <= length(start_nucleo)-idx){
  		i = i + 1
  	}
  	if (i >= 4){
  		write.table(currentstart,file=paste0(x,'Analysis.txt'), append=TRUE, sep='\t', row.names=F, col.names=F, quote=F, eol='\n')
  		# print(c(currentstart, idx,i))
  	}
  }
})

# pdf(file="~/Desktop/plot2.pdf")
# Test <- lapply(files, function(x) {
#   a <- read.table(x)
#   b = a[a[,1]>=-10 & a[,1]<=2020,]
#   c = b[,2]/b[,3]
#   # d = data.frame(x=b[,1],y=c)
#   d = data.frame(x=b[,1],y=c)
#   d$smooth = rollmean(d$y,25,fill='extend')
#   p=ggplot(data=d,aes(x=x,y=smooth))
#   print(p+geom_line()+labs(title=x))
# })
# dev.off()

# print(p+geom_smooth(method = "lm",formula = y ~ poly(x, 8)) +
# coord_cartesian(ylim=c(0,0.1))+labs(title=x))