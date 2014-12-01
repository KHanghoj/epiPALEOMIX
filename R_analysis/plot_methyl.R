require(ggplot2)
require(zoo)
files = list.files(pattern='*Saq*')
pdf(file="~/Desktop/plot2.pdf")
Test <- lapply(files, function(x) {
  a <- read.table(x)
  b = a[a[,1]>=-10 & a[,1]<=2020,]
  c = b[,2]/b[,3]
  # d = data.frame(x=b[,1],y=c)
  d = data.frame(x=b[,1],y=c)
  d$smooth = rollmean(d$y,25,fill='extend')
  p=ggplot(data=d,aes(x=x,y=smooth))
  print(p+geom_line()+labs(title=x))
})
dev.off()

# print(p+geom_smooth(method = "lm",formula = y ~ poly(x, 8)) +
# coord_cartesian(ylim=c(0,0.1))+labs(title=x))