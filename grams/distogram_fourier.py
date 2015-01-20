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
import gzip

_MINMAPQUALI = 25
_MIN_COVERAGE = 1
_OUTLENGTH = int(1e6)  # a million numbers


class Phaso_count():
    def __init__(self):
        self.count = 0


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


def writetofile(output_lst, f_name):
    ''' dfs '''
    with gzip.open('{}.gz'.format(f_name), 'w') as f_output:
        for item in output_lst:
            f_output.write('{}\n'.format(item))


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
        yield (args.chrom, int(args.start), int(args.end))


def filter_by_count(dic):
    return (key for key, value in dic.iteritems() if value >= _MIN_COVERAGE)


def call_output(plus, minus, output_lst, counter_idx):
        plus_good = filter_by_count(plus)
        minus_good = filter_by_count(minus)

        if minus_good and plus_good:
            for plus_pos, minus_pos in product(plus_good, minus_good):
                var = abs(plus_pos-minus_pos)
                output_lst.append(var)
                counter_idx.count += 1


def update(dic, pos):
    try:  # only True if present in dict
        dic[pos] += 1
    except KeyError:
        dic[pos] = 1


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    output_lst = []
    last_tid = -1
    last_pos = -1
    counter_idx = Phaso_count()
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
                call_output(plus, minus, output_lst, counter_idx)
                plus = {}
                minus = {}

            if last_pos < record.pos:  # read not overlapping
                call_output(plus, minus, output_lst, counter_idx)
                plus = {}
                minus = {}

            if record.is_reverse:
                update(minus, record.aend)
            else:
                update(plus, record.pos)
            last_pos = record.aend
            if counter_idx.count > _OUTLENGTH:
                writetofile(output_lst, args.out)
                sys.exit()
        call_output(plus, minus, output_lst, counter_idx)
    writetofile(output_lst, args.out)
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
