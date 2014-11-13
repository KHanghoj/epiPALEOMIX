#!/opt/local/bin/python
'''  Object: return the base composition of the bases in the called nucleosome.
Single-, di- and tetranucleotides '''
from __future__ import print_function
import random
from collections import defaultdict
import argparse
import sys
import pysam
from itertools import product

# Constants
_HALF_NUC = 73
_SINGLENUCLEO = 1
_DINUCLEO = 2
_TETRANUCLEO = 4


def parse_args(argv):
    ''' dfs '''
    parser = argparse.ArgumentParser()
    parser.add_argument('nucleosomepos', help="...")
    parser.add_argument('bam', help="...")
    return parser.parse_args(argv)


def get_dic_fornucl(n_nucl, list_len=147):
    ''' doc '''
    nucleo_dic = defaultdict(int)
    combi = ["".join(map(str, comb)) for comb
             in product('ACGT', repeat=n_nucl)]  # a list of keys
    for i in range(list_len-n_nucl+1):
        nucleo_dic[i] = dict.fromkeys(combi, 0)
    return nucleo_dic


def update_dic(base_list, n_nucl, nucleo_dic):
    ''' doc '''
    str_tempbase = ''.join(base_list)
    for i in range(len(str_tempbase)-n_nucl+1):
        nucleo_dic[i][str_tempbase[i:i+n_nucl]] += 1
    # return nucleo_dic


def get_wind_nucleo(samfile, list_window, chrom, begin_nucl, end_nucl):
    ''' doc '''
    del list_window[:]
    for pileupcolumn in samfile.pileup(chrom, begin_nucl, end_nucl,
                                       truncate=True):

        bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
                 if not x.indel]

        if not bases:
            list_window.append('N')
        else:
            list_window.append(random.choice(bases))


def main(argv):
    ''' dfs '''
    tempbasewin = []
    args = parse_args(argv)
    print(args.nucleosomepos)
    samfile = pysam.Samfile(args.bam, "rb")
    dinucl_update = get_dic_fornucl(_DINUCLEO)
    tetranucl_update = get_dic_fornucl(_TETRANUCLEO)
    singnucl_update = get_dic_fornucl(_SINGLENUCLEO)
    begin_nucl = -1
    end_nucl = -1

    with open(args.nucleosomepos, 'r') as nucleosome_f:
        for line in nucleosome_f.readlines():

            chrom, start, end = (line.rstrip('\n')).split(' ')[:3]

            start = int(start)
            end = int(end)
            dyad = ((end-start)/2) + start
            if ((dyad - _HALF_NUC) == begin_nucl) and ((dyad + _HALF_NUC) ==
                                                       end_nucl):
                continue
            print(chrom, start, end)
            begin_nucl = dyad - _HALF_NUC
            end_nucl = dyad + _HALF_NUC

            get_wind_nucleo(samfile, tempbasewin, chrom, begin_nucl, end_nucl)

            if len(tempbasewin) != 147-1:  # becuase 0-based
                continue

            update_dic(tempbasewin, _SINGLENUCLEO, singnucl_update)
            update_dic(tempbasewin, _DINUCLEO, dinucl_update)
            update_dic(tempbasewin, _TETRANUCLEO, tetranucl_update)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
