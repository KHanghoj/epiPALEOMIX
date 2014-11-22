#!/opt/local/bin/python
'''  Object: To calculate the distograms between mapped reads's start positions
aligning in opposing orientation.
~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --start 16056601 --end 16058000
~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --start 16050500 --end 16050720 --out testdist
16050500 16050720

'''

from __future__ import print_function
from itertools import product
import sys
import pysam
import argparse


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_distogram.txt')
    return parser.parse_args(argv)


def call_distance(l_minus, l_plus, samfile, chrom, f_output):
    ''' docstring '''
    if len(l_minus) > 1 and len(l_plus) > 1:
        # startstrand = True
        # plus = [k for k, v in dic_start.items() if v != startstrand]
        # minus = [k for k, v in dic_start.items() if v == startstrand]
        # create all combinations
        print(l_minus)
        print(l_plus)
        for plus_pos, minus_pos in product(set(l_plus), set(l_minus)):
            print(samfile.getrname(chrom), plus_pos+1, minus_pos+1,
                  abs(plus_pos-minus_pos), file=f_output, sep='\t')
    del l_minus[:]
    del l_plus[:]
    # dic_start.clear()
    #     print(plus, 'plus')
    #     print(minus, 'minus')
    # for _ in sorted(dic_start.keys()):
    #     if len(set(dic_start.values())) > 1:
    #         #   only if data available from both strands
    #         start_pos = min(dic_start)  # start position
    #         startstrand = dic_start[min(dic_start)]  # start position strand
    #         end_pos = max([k for k, v in dic_start.items() if v != startstrand])
    #         length = end_pos-start_pos

    #         print(samfile.getrname(chrom), start_pos+1, end_pos+1,
    #               length, startstrand, (not startstrand), file=f_output, sep='\t')
    #         dic_start.pop(start_pos, None)
    # STILL NEED TO CALCULATE LAST EXAMPLE IN VALOEUV 2011
    # need to pop only smallest value and recalculate


def update_dic(l_minus, l_plus, beginpos, endpos, strand):
    ''' docstring '''
    if strand:  # minus strand
        l_minus.append(endpos)
    else:
        l_plus.append(beginpos)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    f_output = open(args.out, 'w')  # the output file
    chrom = ''
    l_minus = []
    l_plus = []

    for record in samfile.fetch(args.chrom, args.start, args.end):
        strand = record.is_reverse
        if record.tid != chrom:  # new chromosome or first read
            chrom = record.tid
            last_end_pos = record.aend
        if last_end_pos >= record.pos:  # if new read overlaps former
            last_end_pos = record.aend  # assign new end read
            update_dic(l_minus, l_plus, record.pos, last_end_pos, strand)
        else:  # calculate the distance
            call_distance(l_minus, l_plus, samfile, chrom, f_output)
            last_end_pos = record.aend
            update_dic(l_minus, l_plus, record.pos, last_end_pos, strand)

    call_distance(l_minus, l_plus, samfile, chrom, f_output)
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
