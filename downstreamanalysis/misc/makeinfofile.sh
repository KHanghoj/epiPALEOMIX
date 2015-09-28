# takes UCSC as input genes
makePROMGEBO_INFOFILE.R
# ~/data/bedfiles # takes the output of makePROMGEBO_INFOFILE.R
python getcpgs.py
python assignPMcategory.py

# in ~/data/bedfiles
#paste PROM_GEBO_autosom_wochr.INFOFILE <(cut -f 2,3 PROM_GEBO_autosom_wochr.gccontent) <(cut -f 2 PROM_GEBO_autosom_wochr.promcat) > PROM_GEBO_autosom_wochr.INFOFILE.FINAL
#paste PROM_GEBO_autosom_wchr.INFOFILE <(cut -f 2,3 PROM_GEBO_autosom_wochr.gccontent) <(cut -f 2 PROM_GEBO_autosom_wochr.promcat) > PROM_GEBO_autosom_wchr.INFOFILE.FINAL

# add header to these files
#awk 'BEGIN{print "#CHROM\tPROMstart\tPROMend\tGEBOStart\tGEBOend\tSTRAND\tPROMREGION\tGEBOREGION\tOVERALLPROMGC\tOVERALLPROMCpGCOUNT\tPROMCATEGORY"}1' PROM_GEBO_autosom_wochr.INFOFILE.FINAL > PROM_GEBO_autosom_wochr.INFOFILE.FINAL


paste PROM_GEBO_autosom_wochr.INFOFILE <(cut -f 2,3 PROM_GEBO_autosom_wochr.gccontent) <(cut -f 2 PROM_GEBO_autosom_wochr.promcat) |
    awk 'BEGIN{print "#CHROM\tPROMstart\tPROMend\tGEBOStart\tGEBOend\tSTRAND\tPROMREGION\tGEBOREGION\tOVERALLPROMGC\tOVERALLPROMCpGCOUNT\tPROMCATEGORY"}1' > PROM_GEBO_autosom_wochr.INFOFILE.FINAL


paste PROM_GEBO_autosom_wchr.INFOFILE <(cut -f 2,3 PROM_GEBO_autosom_wochr.gccontent) <(cut -f 2 PROM_GEBO_autosom_wochr.promcat) |
    awk 'BEGIN{print "#CHROM\tPROMstart\tPROMend\tGEBOStart\tGEBOend\tSTRAND\tPROMREGION\tGEBOREGION\tOVERALLPROMGC\tOVERALLPROMCpGCOUNT\tPROMCATEGORY"}1' > PROM_GEBO_autosom_wchr.INFOFILE.FINAL
