#!/opt/local/bin/python
'''  Object: return the base composition of the bases in the called nucleosome.
Single-, di- and tetranucleotides '''
from __future__ import print_function

import random
import collections
import argparse

# from collections import defaultdict
# def C(s):
#     dinucleot = 2
#     d = defaultdict(int)
#     # d = {}
#     for i in xrange(len(s)-dinucleot-1):
#         d[s[i:i+dinucleot]] += 1
#     return d

# print(C('ACGTAGCTACGACTCGATGCATGTAGCTAGCTAGCTACGTCAGTAGCTGACTGATCGATCGTCGTAGCTGACTGATCGATCGATCGAATTCCGG'))


    basecomposi = []  # single nucleotide base composition
    cnt = collections.Counter()  # this is for counting bases in basecomposi


# NOTE: Do not include positions that are indels
bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
         if not x.indel]
# NOTE: Don't assume that base 0 is the consensus
# get the consensus nucleotide at each site
# TODO: Handle deletions by checking if bases is empty
# base = random.choice(bases)
s_depth = len(bases)  # get the depth

#
# def parse_args(argv):
#     parser = argparse.ArgumentParser()
#     parser.add_argument('name', help="...")
#     parser.add_argument('--param', help="...")
#     return parser.parse_args(argv)

# def main(argv):
#     args = parse_args(argv)
#     print args.name
#     print args.param

#     return 0

# if __name__ == '__main__':
#     sys.exit(main(sys.argv[1:]))