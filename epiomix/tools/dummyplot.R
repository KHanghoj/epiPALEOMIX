args <- commandArgs(trailingOnly = TRUE)
IN_path = args[1]
outputplot = args[2]
pdf(outputplot)
plot(1, type="n", axes=F, xlab=sprintf("NODATA in %s",IN_path), ylab="")
dev.off()