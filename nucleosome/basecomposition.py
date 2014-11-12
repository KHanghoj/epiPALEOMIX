#!/opt/local/bin/python
'''  Object: return the base composition of the bases in the called nucleosome.
Single-, di- and tetranucleotides '''
from __future__ import print_function

import random
from collections import defaultdict
import argparse
import sys
import pysam

# _NUCLEOSOMERANGE = range(-73, 74)
_HALF_NUC=73


def parse_args(argv):
    ''' dfs '''
    parser = argparse.ArgumentParser()
    parser.add_argument('nucleosomepos', help="...")
    parser.add_argument('bam', help="...")
    return parser.parse_args(argv)


def C(s, length_nuc):
    ''' dingdong '''
    d = defaultdict(int)
    for i in xrange(len(s)-length_nuc-1):
        d[s[i:i+length_nuc]] += 1
    return d


def main(argv):
    ''' dfs '''
    tempbasewin = []
    args = parse_args(argv)
    print(args.nucleosomepos)
    nestbase = []
    samfile = pysam.Samfile(args.bam, "rb")
    with open(args.nucleosomepos, 'r') as f:
        for line in f.readlines():
            chrom, start, end = (line.rstrip('\n')).split(' ')[:3]
            start = int(start)
            end = int(end)
            dyad = ((end-start)/2) + start
            del tempbasewin[:]

            begin_nucl = dyad - _HALF_NUC
            end_nucl = dyad + _HALF_NUC

            for pileupcolumn in samfile.pileup(chrom, begin_nucl, end_nucl,
                                               truncate=True):

                # print(pileupcolumn)
                bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
                         if not x.indel]
                if not bases:
                    tempbasewin.append('N')
                else:
                    # a = random.choice(bases)
                    tempbasewin.append(random.choice(bases))
            # here we need to analyze the bases for the 147: position wize
            nestbase.append(tempbasewin)
        d = defaultdict(int)

        for i in range(10):
            d[i] = {'A': 0, 'C': 0, 'T': 0, 'G' : 0}
        print(d)

        s = 'ACTGACTG'
        s1 = 'ACTGACTG'[::-1]

        for i,h in enumerate(s):
            d[i][h]+=1

        for i,h in enumerate(s1):
            d[i][h]+=1
        print(d)
        # for index, position in enumerate(nestbase-2+1):
            # d[position[index]] += 1
            # print(hmm)
            # print(tempbasewin, len(tempbasewin), sites)

            # for_fetch = (x+dyad for x in _NUCLEOSOMERANGE)
            # print(for_fetch)
            # sys.exit('stop  :)')
        # print(nestbase[1][:])
    return 0






# # print(C('ACGTAGCTACGACTCGATGCATGTAGCTAGCTAGCTACGTCAGTAGCTGACTGATCGATCGTCGTAGCTGACTGATCGATCGATCGAATTCCGG'))

#     # result = list(x.n for x in samfile.pileup())

#     basecomposi = []  # single nucleotide base composition
#     cnt = collections.Counter()  # this is for counting bases in basecomposi


# # NOTE: Do not include positions that are indels
# bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
#          if not x.indel]
# # NOTE: Don't assume that base 0 is the consensus
# # get the consensus nucleotide at each site
# # TODO: Handle deletions by checking if bases is empty
# # base = random.choice(bases)
# s_depth = len(bases)  # get the depth


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
