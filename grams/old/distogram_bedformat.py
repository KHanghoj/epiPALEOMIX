#!/opt/local/bin/python
'''  Object: To calculate the distograms between mapped reads's start positions
aligning in opposing orientation.
~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --start 16056601 --end 16058000
~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --start 16050500 --end 16050720 --out testdist
16050500 16050720

'''

from __future__ import print_function
from itertools import product
from collections import defaultdict
import sys
import pysam
import argparse


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_distogram.txt')
    return parser.parse_args(argv)


# def call_distance(l_minus, l_plus, samfile, chrom, f_output):
def call_distance(l_minus, l_plus, output_dic):
    ''' docstring '''
    if len(l_minus) > 0 and len(l_plus) > 0:
        for plus_pos, minus_pos in product(set(l_plus), set(l_minus)):
            var = abs(plus_pos-minus_pos)
            output_dic[var] += 1  # this creates much smaller output files.
            # print(samfile.getrname(chrom), plus_pos+1, minus_pos+1,
                  # abs(plus_pos-minus_pos), file=f_output, sep='\t')
    del l_minus[:]
    del l_plus[:]


def update_dic(l_minus, l_plus, beginpos, endpos, strand):
    ''' docstring '''
    if strand:  # minus strand
        l_minus.append(endpos)
    else:
        l_plus.append(beginpos)


def writetofile(output_dic, f_output):
    for key, value in output_dic.iteritems():
        f_output.write('{}\t{}\n'.format(key, value))


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    f_output = open(args.out, 'w')  # the output file
    bedfile = args.bed
    chrom = None
    output_dic = defaultdict(int)
    # make a class for plus strand and minus strand with
    # positions and pileup count.
    l_minus = []
    l_plus = []
    with open(bedfile, 'r') as bedfile_f:
        for line in bedfile_f.readlines():
            input_line = (line.rstrip('\n')).split('\t')[:3]

            # it is the users responsibility to input bed format
            # identical to BAM format.
            try:
                chrom = input_line.pop(0).replace('chr', '')
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
            except ValueError:
                chrom = args.chrom
                start = args.start
                end = args.end

            for record in samfile.fetch(chrom, start, end):
                if record.tid != chrom:  # new chromosome or first read
                    chrom = record.tid
                    last_end_pos = record.aend
                if last_end_pos >= record.pos:  # if new read overlaps former
                    last_end_pos = record.aend  # assign new end read
                    update_dic(l_minus, l_plus, record.pos,
                               last_end_pos, record.is_reverse)
                else:  # calculate the distance and restart
                    call_distance(l_minus, l_plus, output_dic)
                    last_end_pos = record.aend
                    update_dic(l_minus, l_plus, record.pos,
                               last_end_pos, record.is_reverse)

    call_distance(l_minus, l_plus, output_dic)
    writetofile(output_dic, f_output)
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
