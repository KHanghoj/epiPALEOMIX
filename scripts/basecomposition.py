#!/opt/local/bin/python
'''  Object: return the base composition of the bases in the called nucleosome.
Single-, di- and tetranucleotides
### note: we need to find a solution of nucleosome mapping with wide peaks.
### need to change the code to fit only one dyad
python ~/install/python_scripts/basecomposition.py ../subsamples/nucleomap/subsample_0_80_Saqqaq ../Saqqaq.hg19.flt.sort.rmdup.realign.md.bam
'''
from __future__ import print_function
from random import choice
from collections import defaultdict
import argparse
import sys
import pysam
from itertools import product

# Constants
# _HALF_NUC = 73  #
_HALF_NUC = 120  # a bit larger because of genome aresearch paper.
_SINGLENUCLEO = 1
_DINUCLEO = 2
_TETRANUCLEO = 4
_SCORE = 2


def parse_args(argv):
    ''' dfs '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="special bedformat (the output from\
                        the nucleomap.py script")
    parser.add_argument('--out', help='...', default='out_basecomposition.txt')
    return parser.parse_args(argv)


def get_dic_fornucl(n_nucl, list_len=147):
    ''' doc '''
    nucleo_dic = defaultdict(int)
    combi = ["".join(map(str, comb)) for comb
             in product('ACGT', repeat=n_nucl)]  # a list of keys
    for i in range(0, list_len-n_nucl+1, n_nucl):
        nucleo_dic[i] = dict.fromkeys(combi, 0)
    return nucleo_dic


def update_dic(base_list, n_nucl, nucleo_dic):
    ''' doc '''
    str_tempbase = ''.join(base_list)
    for i in range(0, len(base_list)-n_nucl+1, n_nucl):
    # extra minus 1 is because the len() is not zero based
        try:
            nucleo_dic[i][str_tempbase[i:i+n_nucl]] += 1
        except KeyError:
            continue
    # return nucleo_dic


def get_wind_nucleo(samfile, tempbasewin, chrom, begin_nucl, end_nucl):
    ''' doc '''
    del tempbasewin[:]
    for pileupcolumn in samfile.pileup(chrom, begin_nucl, end_nucl,  # ):
                                       truncate=True):
        bases = [x.alignment.seq[x.query_position] for x in pileupcolumn.pileups
                 if not x.indel]
        if not bases:
            tempbasewin.append('N')
        else:
            tempbasewin.append(choice(bases))


def writetofile(dic, f_name, args):
    ''' dfs '''
    full_f_name = '{}{}'.format(args.out, f_name)
    f_output = open(full_f_name, 'w')
    word_combinations = [x for x in sorted(dic[0].keys())]
    f_output.write('\t'+'{}'.format('\t'.join(word_combinations))+'\n')
    for i in sorted(dic.keys()):
        f_output.write(repr(i))
        for value in sorted(dic[i]):  # sort word combinations. A,C,G,T
            # f_output.write('\t{}:{}'.format(repr(value), repr(dic[i][value])))
            f_output.write('\t{}'.format(repr(dic[i][value])))
        f_output.write('\n')

    f_output.close()


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')
                chrom = input_line.pop(0)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                depth = int(input_line.pop(0))
                score = float(input_line.pop(0))
                yield (chrom, start, end, depth, score)
    else:
        yield (args.chrom, args.start, args.end)


def main(argv):
    ''' dfs '''
    tempbasewin = []
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")

    begin_nucl = -1
    end_nucl = -1
    last_chrom = -1
    last_dyad = -1
    windowsize = len(range(0, _HALF_NUC*2+1))

    dinucl_update = get_dic_fornucl(_DINUCLEO, windowsize)
    tetranucl_update = get_dic_fornucl(_TETRANUCLEO, windowsize)
    singnucl_update = get_dic_fornucl(_SINGLENUCLEO, windowsize)
    # with open(args.nucleosomepos, 'r') as nucleosome_f:
    #     for line in nucleosome_f.readlines():
    for chrom, start, end, depth, score in read_bed(args):
            if score < _SCORE:  # arbitrary value
                continue
            if abs(start-end) > 0:  # remove wide nuclesome calls
                continue
            dyad = start
            if dyad == last_dyad and chrom == last_chrom:
                continue  # will not use the same position twice
            last_chrom = chrom
            last_dyad = dyad
            begin_nucl = dyad - _HALF_NUC
            end_nucl = dyad + _HALF_NUC + 1
            get_wind_nucleo(samfile, tempbasewin, chrom, begin_nucl, end_nucl)
            if len(tempbasewin) != windowsize:  # becuase 0-based
                continue
            update_dic(tempbasewin, _SINGLENUCLEO, singnucl_update)
            update_dic(tempbasewin, _DINUCLEO, dinucl_update)
            update_dic(tempbasewin, _TETRANUCLEO, tetranucl_update)
    writetofile(singnucl_update, 'single', args)
    writetofile(dinucl_update, 'dinucleo', args)
    writetofile(tetranucl_update, 'tetra', args)
    samfile.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))