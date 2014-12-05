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
_MIN_COVERAGE = 5


class Phaso():
    def __init__(self, pos):
        self.position = int(pos)
        self.count = 0


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', help="...")
    parser.add_argument('bam', help="...")
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
    while starts[-1].position - starts[0].position > max_lst_range:
        oldest = starts.pop(0)
        for current in starts:
            # Have we reached positions outside
            # the maximum range we are interested in?
            if current.position - oldest.position >= max_size:
                break

            count = oldest.count
            length = current.position - oldest.position
            if count >= _MIN_COVERAGE:
                output_dic[length] += 1
                # print count, length


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
    starts = [Phaso(0)]  # initialize the positions and counts
    ends = [Phaso(0)]
    last_tid = -1

    for chrom, start, end in read_bed(args):
        starts = [Phaso(0)]
        ends = [Phaso(0)]

        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            if last_tid != record.tid:
                last_tid = record.tid
                call_output(starts, output_dic, max_lst_range=0)
                call_output(ends, output_dic, max_lst_range=0)
                starts = [Phaso(0)]
                ends = [Phaso(0)]

            if record.is_reverse:
                pos, lst = record.aend, ends
            else:
                pos, lst = record.pos, starts

            if lst[-1].position != pos:
                lst.append(Phaso(pos))
                call_output(lst, output_dic)
            lst[-1].count += 1
        call_output(ends, output_dic, max_lst_range=0)
        call_output(starts, output_dic, max_lst_range=0)

    writetofile(output_dic, args.out)
    samfile.close()
    return 0

    # it = read_bed(args)
    # while True:  # while true needs to have a break within
                   # or if in a function, it needs a break
                   # outside the function if in a loop
                   # like below
    #     try:
    #         chrom, start, end = it.next()

    #         ....
    #     except StopIteration:
    #         break

    # def infinite():
    #     x = 0
    #     while True:
    #         yield x
    #         x += 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
