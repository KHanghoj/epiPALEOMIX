PROMwo <- read.table('/Users/kehanghoej/newtemp/PROM_autosom_wochr.bed')
PROMw <- read.table('/Users/kehanghoej/newtemp/PROM_autosom_wchr.bed')
GEBOwo <- read.table('/Users/kehanghoej/newtemp/GEBO_autosom_wochr.bed')
GEBOw <- read.table('/Users/kehanghoej/newtemp/GEBO_autosom_wchr.bed')
promcatwo <- read.table('/Users/kehanghoej/newtemp/fastaseqs/PROM_autosom_wochr.promcat',h=T)
promcatw <- promcatwo
promcatw$name <- with(promcatw, paste0('chr', name))



printfunc <- function(prom, GEBO, dfcat, ending){
    write.table(prom[with(dfcat, PromGroup=='HIGH'),],
                file=sprintf('/Users/kehanghoej/newtemp/PROM_autosom_HIGH_%s.bed',ending),
                row.names=F,col.names=F,quote=F,sep='\t')
    write.table(prom[with(dfcat, PromGroup=='INTERMEDIATE'),],
                file=sprintf('/Users/kehanghoej/newtemp/PROM_autosom_INTERMEDIATE_%s.bed', ending),
                row.names=F,col.names=F,quote=F,sep='\t')
    write.table(prom[with(dfcat, PromGroup=='LOW'),],
                file=sprintf('/Users/kehanghoej/newtemp/PROM_autosom_LOW_%s.bed',ending),
                row.names=F,col.names=F,quote=F,sep='\t')

    write.table(GEBO[with(dfcat, PromGroup=='HIGH'),],
                file=sprintf('/Users/kehanghoej/newtemp/GEBO_autosom_HIGH_%s.bed',ending),
                row.names=F,col.names=F,quote=F,sep='\t')
    write.table(GEBO[with(dfcat, PromGroup=='INTERMEDIATE'),],
                file=sprintf('/Users/kehanghoej/newtemp/GEBO_autosom_INTERMEDIATE_%s.bed', ending),
                row.names=F,col.names=F,quote=F,sep='\t')
    write.table(GEBO[with(dfcat, PromGroup=='LOW'),],
                file=sprintf('/Users/kehanghoej/newtemp/GEBO_autosom_LOW_%s.bed',ending),
                row.names=F,col.names=F,quote=F,sep='\t')

}

printfunc(PROMwo, GEBOwo, promcatwo, 'wochr')
printfunc(PROMw, GEBOw,promcatw, 'wchr')
                                        # promcatw not really needed.

