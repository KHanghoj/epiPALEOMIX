from __future__ import print_function
from epipaleomix.tools.commonutils import Cache
import os.path, re
#files = ['PROM_minusstrand_autosom_wochr.fa', 'PROM_plusstrand_autosom_wochr.fa']
files = ['PROM_autosom_wochr.fa']
WINDOWLENGTH = 500
HIGHCPG = 0.75
MIDCPG = 0.48
HIGHCGs = 0.55
FMT = '{}\t{}\n'.format
referencepath='/home/krishang/data/reference_human/hs.build37.1.fa'

def unpack(chrom, start, end, *rest):
    return str(chrom), int(start), int(end)

def CpGyield(seq):
    l = len(seq)
    for idx in xrange(l-WINDOWLENGTH+1):
        dat = seq[idx:(idx+WINDOWLENGTH)]
        cgs = (dat.count('C')+dat.count('G'))/float(WINDOWLENGTH)
        exp = (cgs/2)*(cgs/2)
        obs = dat.count('CG')/float(WINDOWLENGTH)
        yield cgs, obs/exp

def getcpgs(seq):
    flag = 'LOW'
    for cgs, cpg in CpGyield(seq):
        if cgs >= HIGHCGs and cpg >= HIGHCPG:
            return 'HIGH'
        elif cpg >= MIDCPG:
            flag = 'INTERMEDIATE'
    return flag

for fin in files:
    fasta = Cache(referencepath)
    with open(fin,'r') as infile:
        fout, _ = os.path.splitext(fin)
        with open(fout+'.promcat', 'w') as outfile:
            for line in infile:
                chrom, start, end = unpack(*re.split(r'\s', line.rstrip()))
                regionname = '{}_{}_{}'.format(chrom,start, end)
                # start-1 because i wanna include the actual bases at position start.
                # bedfiles are always 1-based. fasta.fetch_string is 0based
                startzero = start-1
                outfile.write(FMT(regionname, getcpgs(fasta.fetch_string(chrom, startzero, end-startzero))))
