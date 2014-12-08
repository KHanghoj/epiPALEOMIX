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


# _MAX_SIZE = 147
_MINMAPQUALI = 25
_MIN_COVERAGE = 3


class Phaso():
    def __init__(self, pos):
        self.position = int(pos)
        self.count = 0
        # self.dic = {int(pos): 0}


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
    f_output = open(f_name, 'w')
    for key in sorted(output_dic):
        f_output.write('{}\t{}\n'.format(key, output_dic[key]))
    f_output.close()


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


def filter_by_count(lst):
    return [x.position for x in lst if x.count >= _MIN_COVERAGE]


def call_output(plus, minus, output_dic):
        plus_good = filter_by_count(plus)
        minus_good = filter_by_count(minus)

        if minus_good and plus_good:
            for plus_pos, minus_pos in product(plus_good, minus_good):
                var = abs(plus_pos-minus_pos)
                output_dic[var] += 1  # this creates much smaller output files.


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    output_dic = defaultdict(int)
    plus = [Phaso(0)]  # initialize the positions and counts
    minus = [Phaso(0)]
    last_tid = -1
    last_pos = -1

    for chrom, start, end in read_bed(args):
        plus = [Phaso(0)]
        minus = [Phaso(0)]
        last_pos = -1
        last_tid = -1

        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            if last_tid != record.tid:
                last_tid = record.tid
                last_pos = record.aend
                call_output(plus, minus, output_dic)
                plus = [Phaso(0)]
                minus = [Phaso(0)]

            if last_pos < record.pos:  # read not overlapping
                call_output(plus, minus, output_dic)
                plus = [Phaso(0)]
                minus = [Phaso(0)]

            THINK ABOUT MAKING A DICT INSTEAD OF THE LIST OF POS. SO MUCH FASTER FOR MINUS STRAND THAT IS UNSORTED.

            if record.is_reverse:
                # if minus[-1].position != record.aend:
                    # minus.append(Phaso(record.aend))
                # minus[-1].count += 1
                try:
                    idx = [x.position for x in minus].index(record.aend)
                    minus[idx].count += 1
                except ValueError:
                    minus.append(Phaso(record.aend))
                    minus[-1].count += 1
                # pos, lst = record.aend, minus
            else:
                if plus[-1].position != record.pos:
                    plus.append(Phaso(record.pos))
                plus[-1].count += 1
                # pos, lst = record.pos, plus

            # if lst[-1].position != pos:
                # lst.append(Phaso(pos))
            # lst[-1].count += 1

            # if lst[-1].position != pos:
            #     lst.append(Phaso(pos))
            # lst[-1].count += 1

            # if record.is_reverse:
            #     try:

            #     except ValueError:
            #         print('adsfasdfasdf')
            #         print(lst[:].count)
            # print()

            # call_output(lst, output_dic)
            last_pos = record.aend
        call_output(plus, minus, output_dic)
        # call_output(ends, output_dic, max_lst_range=0)
        # call_output(starts, output_dic, max_lst_range=0)
    call_output(plus, minus, output_dic)
    writetofile(output_dic, args.out)
    samfile.close()
            # if last_end_pos >= record.pos:  # if new read overlaps former
    #             last_end_pos = record.aend  # assign new end read
    #             update_dic(l_minus, l_plus, record.pos,
    #                        last_end_pos, record.is_reverse)
    #         else:  # calculate the distance and restart
    #             call_distance(l_minus, l_plus, output_dic)
    #             last_end_pos = record.aend
    #             update_dic(l_minus, l_plus, record.pos,
    #                        last_end_pos, record.is_reverse)

    # call_distance(l_minus, l_plus, output_dic)
    # writetofile(output_dic, f_output)
    # f_output.close()
    # samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
