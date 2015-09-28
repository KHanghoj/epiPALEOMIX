from __future__ import print_function
import pysam, re, sys
fastpath = '/home/krishang/data/reference_human/hs.build37.1.fa'
#fastpath = 'chrome.fa'
CGpattern = re.compile('CG')
FMT = '{}\t{}\t{}\n'.format

def getindex(seq, pat=CGpattern):
    for m in pat.finditer(seq):
        yield m.start()

def getchroms(path):
    with open(path+'.fai') as f_chrom:
        return [line.split('\t', 1)[0] for line in f_chrom]


chroms = getchroms(fastpath)
fasta = pysam.Fastafile(fastpath)
for chrom in chroms:
    print(chrom, file=sys.stderr, end='\t') 
    size = fasta.get_reference_length(chrom)
    jump = int(1e7)
    start = 1
    while start < size:
        end = start+jump
        for s in getindex(fasta.fetch(chrom,start-1,end)):
            beg = start+s # this is 1-based for BED format
            sys.stdout.write(FMT(chrom, beg, beg+1))
        start = end
sys.stdout.close()
