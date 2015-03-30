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


_MINMAPQUALI = 25
_MIN_COVERAGE = 3


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_distogram.txt')
    return parser.parse_args(argv)


def writetofile(output_dic, f_name):
    ''' dfs '''
    with open(f_name, 'w') as f_output:
        for key in sorted(output_dic):
            f_output.write('{}\t{}\n'.format(key, output_dic[key]))


# def writenotabulate(notab, out):
#     import gzip
#     outf = 'out_notab{}.gz'.format(out[3:])
#     tab = iter(notab)
#     with gzip.open(outf, 'wb') as f:
#         [f.write('{}\n'.format(val)) for val in tab]
#         # for val in tab:
#             # f.write('{}\n'.format(val))


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')[:3]
                chrom = input_line.pop(0).replace('chr', '')
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        yield (args.chrom, args.start, args.end)


def filter_by_count(dic):
    return (key for key, value in dic.iteritems() if value >= _MIN_COVERAGE)


def call_output(plus, minus, output_dic):
        plus_good = filter_by_count(plus)
        minus_good = filter_by_count(minus)

        if minus_good and plus_good:
            for plus_pos, minus_pos in product(plus_good, minus_good):
                var = abs(plus_pos-minus_pos)
                output_dic[var] += 1  # this creates much smaller output files.


def update(dic, pos):
    try:  # only True if present in dict
        dic[pos] += 1
    except KeyError:
        dic[pos] = 1


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    output_dic = defaultdict(int)
    plus = {}
    minus = {}
    last_tid = -1
    last_pos = -1
    for chrom, start, end in read_bed(args):
        plus = {}
        minus = {}
        last_pos = -1
        last_tid = -1

        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            if last_tid != record.tid:
                last_tid = record.tid
                last_pos = record.aend
                call_output(plus, minus, output_dic)
                plus = {}
                minus = {}

            if last_pos < record.pos:  # read not overlapping
                call_output(plus, minus, output_dic)
                plus = {}
                minus = {}

            if record.is_reverse:
                update(minus, record.aend)
            else:
                update(plus, record.pos)
            last_pos = record.aend
        call_output(plus, minus, output_dic)
    writetofile(output_dic, args.out)
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
