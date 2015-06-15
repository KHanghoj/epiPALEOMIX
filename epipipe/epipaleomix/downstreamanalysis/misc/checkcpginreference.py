from __future__ import print_function
from epipaleomix.tools.commonutils import Cache
import re

def unpack(c,s,e):
    return str(c), int(s), int(e)

referencepath='/home/krishang/data/reference_human/hs.build37.1.fa'
fasta = Cache(referencepath)
with open('allcoordinates.bed', 'r') as f_in:
    with open('allcoordinatescpgchecked.bed', 'w') as f_out:
        for line in f_in:
            chrom, start, end = unpack(*re.split(r'\s', line.rstrip()))
            startzero = start-1  # bed files are 1-based. my Cache is 0-based
            if fasta.fetch_string(chrom, startzero, end-startzero) in 'CG':
                f_out.write(line)
