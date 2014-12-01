require(ggplot2)
require(zoo)
files = list.files(pattern='*Saq*')
pdf(file="plot.pdf")
Test <- lapply(files, function(x) {
  a <- read.table(x)
  b = a[a[,1]>=10 & a[,1]<=1020,]
  # d = data.frame(x=b[,1],y=c)
  d = data.frame(x=b[,1],y=b[,2])
  d$smooth = rollmean(d$y,4,fill='extend')
  p=ggplot(data=d,aes(x=x,y=smooth))
  print(p+geom_bar(stat='identity', binwidth=1)+labs(title=x))
})
dev.off()