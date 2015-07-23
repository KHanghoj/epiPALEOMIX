from __future__ import print_function
import gzip
l = set()
with open('../HK_genes.txt','r') as fcheck:
    for line in fcheck:
        l.add(line.rstrip().split('\t')[1])
        
with open('Rsscorefortest_withoutHK.txt','w') as fout:
    with gzip.open('../Rsscores.txt.gz','rb') as fin:
        for line in fin:
            if not line.rstrip().split('\t')[24] in l:
                fout.write(line)
