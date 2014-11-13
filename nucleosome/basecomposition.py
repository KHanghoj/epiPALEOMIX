#!/opt/local/bin/python
'''  Object: return the base composition of the bases in the called nucleosome.
Single-, di- and tetranucleotides '''
from __future__ import print_function

import random
from collections import defaultdict
import argparse
import sys
import pysam
import itertools


# _NUCLEOSOMERANGE = range(-73, 74)
_HALF_NUC=73


def parse_args(argv):
    ''' dfs '''
    parser = argparse.ArgumentParser()
    parser.add_argument('nucleosomepos', help="...")
    parser.add_argument('bam', help="...")
    return parser.parse_args(argv)


# def C(s, length_nuc):
#     ''' dingdong '''
#     d = defaultdict(int)
#     for i in xrange(len(s)-length_nuc-1):
#         d[s[i:i+length_nuc]] += 1
#     return d


def get_dic_fornucl(n_nucl, list_len=147):
    ''' doc '''
    d = defaultdict(int)
    combi = ["".join(map(str, comb)) for comb
             in itertools.product('ACGT', repeat=n_nucl)]  # a list of keys
    for i in range(list_len-n_nucl+1):
        d[i] = dict.fromkeys(combi, 0)
    return d


def update_dic(base_list, n_nucl, d):
    ''' doc '''
    str_tempbase = ''.join(base_list)
    for i in range(len(str_tempbase)-n_nucl+1):
        d[i][str_tempbase[i:i+n_nucl]] += 1
    return d


def main(argv):
    ''' dfs '''
    tempbasewin = []
    args = parse_args(argv)
    print(args.nucleosomepos)
    samfile = pysam.Samfile(args.bam, "rb")
    singlenucleotides = 1
    dinucletides = 2
    tetranucleotides = 4
    dinucl_update = get_dic_fornucl(dinucletides)
    tetranucl_update = get_dic_fornucl(tetranucleotides)
    singnucl_update = get_dic_fornucl(singlenucleotides)
    begin_nucl = -1
    end_nucl = -1
    sites = 0
    with open(args.nucleosomepos, 'r') as f:
        for line in f.readlines():
            chrom, start, end = (line.rstrip('\n')).split(' ')[:3]

            start = int(start)
            end = int(end)
            dyad = ((end-start)/2) + start
            del tempbasewin[:]
            if ((dyad - _HALF_NUC) == begin_nucl) and ((dyad + _HALF_NUC) ==
                                                       end_nucl):
                continue
            print(chrom, start, end)
            begin_nucl = dyad - _HALF_NUC
            end_nucl = dyad + _HALF_NUC

            for pileupcolumn in samfile.pileup(chrom, begin_nucl, end_nucl,
                                               truncate=True):

                bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
                         if not x.indel]
                if len(set(bases))>1: print (bases) # we have some variants in the data set.
                if not bases:
                    tempbasewin.append('N')
                else:
                    # a = random.choice(bases)
                    tempbasewin.append(random.choice(bases))
            print(len(tempbasewin))
            if len(tempbasewin) != 147-1:  # becuase 0-based
                continue
            print(''.join(tempbasewin))
            print(singnucl_update[0])
            print()
            singnucl_update = update_dic(tempbasewin, singlenucleotides,
                                         singnucl_update)
            dinucl_update = update_dic(tempbasewin, dinucletides, dinucl_update)
            tetranucl_update = update_dic(tempbasewin, tetranucleotides,
                                          tetranucl_update)
            sites += 1

            print(sites)
            print(singnucl_update[0])
            # exit('dsfdsfsd')

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
