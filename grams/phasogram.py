#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict

_MAX_SIZE = 3000
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
    parser.add_argument('--out', help='...', default='out_phasogram.txt')
    return parser.parse_args(argv)


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


def call_output(starts, output_dic,
                max_lst_range=_MAX_SIZE, max_size=_MAX_SIZE):
    if starts:
        while max(starts) - min(starts) > max_lst_range:
            old_pos = min(starts)
            old_count = starts.pop(old_pos, None)
            for current in sorted(starts.iterkeys()):
                length = current - old_pos
                if length >= max_size:
                    break
                    ## do no know if break or contiune
                if old_count >= _MIN_COVERAGE:
                    output_dic[length] += 1


def writetofile(output_dic, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    for key, value in output_dic.iteritems():
        f_output.write('{}\t{}\n'.format(key, value))
    f_output.close()





def main(argv):
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    output_dic = defaultdict(int)
    starts = {}  # initialize the positions and counts
    ends = {}
    last_tid = -1

    for chrom, start, end in read_bed(args):
        starts = {}
        ends = {}

        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            if last_tid != record.tid:
                last_tid = record.tid
                call_output(starts, output_dic, max_lst_range=0)
                call_output(ends, output_dic, max_lst_range=0)
                starts = {}
                ends = {}
            if record.is_reverse:
                pos, present_dic = record.aend, ends
            else:
                pos, present_dic = record.pos, starts

            if present_dic.get(pos, 0):  # only True if present pos in dict
                present_dic[pos] += 1
            else:
                present_dic[pos] = 1
                call_output(present_dic, output_dic)

        call_output(ends, output_dic, max_lst_range=0)
        call_output(starts, output_dic, max_lst_range=0)

    writetofile(output_dic, args.out)
    samfile.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
