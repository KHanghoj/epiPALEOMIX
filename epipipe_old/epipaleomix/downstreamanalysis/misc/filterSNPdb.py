from __future__ import print_function
import gzip, re, sys
from epipaleomix.tools.commonutils import Cache

def unpackSNP(chrom, start, end, snprs, strand, base, derived):
    return {'c':str(chrom), 's':int(start), 'outpos': int(end), 
            'strand':str(strand), 'refbas':str(base), 'derived':str(derived)}

FMT='{c}\t{outpos}\t{strand}\t{refbas}\t{derived}\n'.format
#outname, infilepath = sys.argv[1:]
infilepath, outname = sys.argv[1:]
fastapath = '/home/krishang/data/reference_human/hs.build37.1.fa'
fasta = Cache(fastapath, seq_len=1000)
CG = 'CG'
print(infilepath, outname , file=sys.stderr)
with gzip.open(outname, 'wb') as fout:
    with gzip.open(infilepath, 'rb') as inf:
        ## the start position in the snp database is the base before the snp
        ## so indirectly 0-based. It there fits my 0-based index
        ## no header to remove

## this is for C to T only:
        # for line in inf:
        #     d = unpackSNP(*re.split(r'\s+', line.rstrip()))
        #     # re.match return None if no hit as the beginning
        #     if d['strand'] == '+' and d['refbas'] == 'C' and re.match('C/T', d['derived']):
        #         if fasta.fetch_string(d['c'], d['s'], 2) == CG:
        #             fout.write(FMT(**d))
        #     elif d['strand'] == '-' and d['refbas'] == 'G' and re.match('C/T', d['derived']):
        #         if fasta.fetch_string(d['c'], d['s']-1, 2) == CG:
        #             d['outpos'] -= 1
        #             fout.write(FMT(**d))
## this is for C to N still in CpG context:
        for line in inf:
            d = unpackSNP(*re.split(r'\s+', line.rstrip()))
            # re.match return None if no hit as the beginning
            if d['strand'] == '+' and d['refbas'] == 'C' and fasta.fetch_string(d['c'], d['s'], 2) == CG:
                fout.write(FMT(**d))
            elif d['strand'] == '-' and d['refbas'] == 'G' and fasta.fetch_string(d['c'], d['s']-1, 2) == CG:
                d['outpos'] -= 1
                fout.write(FMT(**d))
