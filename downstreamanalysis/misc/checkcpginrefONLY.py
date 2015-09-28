from __future__ import print_function
from epipaleomix.tools.commonutils import Cache
import fileinput
import re, sys

def unpack(c,s,*rest):
    return str(c), int(s)

## remember to sort files before checking
## sort -k 1 -k 2n Hrce_1.bed| head
## sort -k 1 -k 2n Hrce_1.bed|python  ../../pythonscripts/checkcpginref.py > Hrce_1_checked.bed

referencepath = '/home/krishang/data/reference_human/hs.build37.1.fa'
# referencepath='/home/krishang/data/reference_human/hs.build37.chr1.fa'
fasta = Cache(referencepath, seq_len=1000)
checked = 0
for line in fileinput.input():
    if checked % 20000 == 0:
        print(checked, file=sys.stderr,end='\t')
    chrom, start = unpack(*re.split(r'\s', line.rstrip()))
    startzero = start-1
    dinucl = fasta.fetch_string(chrom, startzero, 2)
    if dinucl and dinucl in 'CG':
        sys.stdout.write(line)
    checked +=1
