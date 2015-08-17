# !/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
'''
from __future__ import print_function
import sys
import pysam
import argparse
import re
import gzip
from os.path import exists, splitext
from shutil import move
from collections import defaultdict, namedtuple
from epipaleomix.tools.commonutils import \
    read_bed, \
    Cache
CONV_REV = {'N': 'N',
            'C': 'G',
            'T': 'A',
            'G': 'C',
            'A': 'T'}
CONV_FORW = {'N': 'N',
             'C': 'C',
             'T': 'T',
             'G': 'G',
             'A': 'A'}
revstrand = {True: 'MinusStrand', False:'PlusStrand'}

## _fasta = Cache('/home/krishang/data/reference_human/hs.build37.1.fa')
_fasta = Cache('chrome.fa')
ReadBases=2
chrom='22'
pat = re.compile('CG')

def _getindexes(bases_str):
    ''' returns the 0-based indeces of fasta read'''
    return (m.start() for m in pat.finditer(bases_str))

def _prep(curr_pos, skip=0):
    ## skip is only used in right side functions
    # skip is 1-based always
    fast_string = _fasta.fetch_string(chrom, curr_pos,
                                           ReadBases-skip)
    for fast_idx in _getindexes(fast_string):
        yield fast_idx, ReadBases - fast_idx

def _rightpart(record, skip=0):
    curr_pos = record.aend-ReadBases
    cigar_op, cigar_len = record.cigar[-1]
    bases = record.seq[-ReadBases:]
    for fast_idx, inverse_idx in _prep(curr_pos, skip):
        if (cigar_op == 0 and cigar_len >= inverse_idx):
            #return (bases[fast_idx:(fast_idx+2)], record.get_tag('RG'))
            return bases[fast_idx:(fast_idx+2)]

def _leftpart(record, skip=0):
    curr_pos = record.pos
    cigar_op, cigar_len = record.cigar[0]
    bases = record.seq[:ReadBases]
    for fast_idx, _ in _prep(curr_pos):
        if (fast_idx >= skip and cigar_op == 0 and
            cigar_len >= fast_idx+2):
            #return (bases[fast_idx:(fast_idx+2)], record.get_tag('RG'))
            return bases[fast_idx:(fast_idx+2)]
f=str(sys.argv[1])
##records= pysam.AlignmentFile('AltaiNea_chrom22snip.bam')
##records= pysam.AlignmentFile('saqqaq_chrom22.bam').fetch('22', 20000000,25000000)
##records= pysam.AlignmentFile('Losch_chrom22_snip.bam').fetch('22', 20000000,25000000)
records= pysam.AlignmentFile(f).fetch('22', 20000000,25000000)
LIB='Loschbour'
LIB=['SL3003','SL3004']
LIB = sys.argv[2:]
print(f, LIB, file=sys.stderr)
from collections import Counter


forward_five, forward_three= [],[]
reverse_five, reverse_three = [], []

for record in records:
    if record.get_tag('RG') not in LIB:
        if record.is_reverse:
            if _fasta.fetch_string('22', record.pos,2 ) == 'CG':
                reverse_three.append(_leftpart(record))
            if _fasta.fetch_string('22', record.aend-2,2 ) == 'CG':
                reverse_five.append(_rightpart(record))
            
        else:
            if _fasta.fetch_string('22', record.pos,2 ) == 'CG':
                forward_five.append(_leftpart(record))
            if _fasta.fetch_string('22', record.aend-2,2 ) == 'CG':
                forward_three.append(_rightpart(record))

# print('Forward strand; Five prime: ', Counter(forward_five))
# print('Forward strand; Three prime: ', Counter(forward_three))
# print('Reverse strand; Five prime: ', Counter(reverse_five))
# print('Reverse strand; Three prime: ', Counter(reverse_three))

def printout(lst, name):
    for di, coun in Counter(lst).iteritems():
        print(di, coun, name,sep='\t',file=sys.stdout)


printout(forward_five, 'forward_five')
printout(forward_three, 'forward_three')
printout(reverse_five, 'reverse_five')
printout(reverse_three, 'reverse_three')

# for di, coun in Counter(forward_five).iteritems():
#     print(di, coun, 'forward_five',sep='\t',file=sys.stdout)
