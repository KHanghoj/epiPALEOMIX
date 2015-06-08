saq <- read.table('Saqqaq_MethylMap_PROMHIGH_bedcoord.txt.gz')
ust <- read.table('UstIshim_MethylMap_PROMHIGH_bedcoord.txt.gz')
saq$V10 <- with(saq, sprintf('%s_%s_%s', V1,V2,V3))
ust$V10 <- with(ust, sprintf('%s_%s_%s', V1,V2,V3))
df <- merge(saq,ust, by='V10', suffixes = c(".saq",".ust"))
df$delta <- with(df,V5.saq-V5.ust)

alt <- read.table('/Users/kehanghoej/Desktop/AltaiNeanderthal_MethylMap_PROMHIGH_bedcoord.txt.gz')
den <- read.table('/Users/kehanghoej/Desktop/Denisova_MethylMap_PROMHIGH_bedcoord.txt.gz')
alt$V10 <- with(alt, sprintf('%s_%s_%s', V1,V2,V3))
den$V10 <- with(den, sprintf('%s_%s_%s', V1,V2,V3))
df <- merge(alt,den, by='V10', suffixes = c(".alt",".den"))
df$delta <- with(df,V5.alt-V5.den)

summary(lm(V9.alt~V9.den*delta  , data=df))




load('complete450kdatasetwithcoord_sorted.Rdata')  #CALLED metdata  samples from column 2 to 57 incl.
a <- read.table('Saqqaq_MethylMap_ALLPROM_hashremoved_bedcoord.txt')
df <- (a[a$V8>60,])

getrows <- function(rowidx){
    row <- unlist(a[rowidx,1:3])
    md <- metdata[metdata$chr==as.character(row[1]),]
    md <- md[md$genomicpos>=row[2],]
    c(row, colMeans(md[md$genomicpos<=row[3],2:57]))
}
dfs <- do.call(rbind, parallel::mclapply(1:nrow(a),getrows,mc.cores=3))

dim(dfs)
dim(a)

cbind()

1  860029  861529

