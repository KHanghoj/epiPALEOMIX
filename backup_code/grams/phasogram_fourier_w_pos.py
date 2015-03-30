#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
import gzip

_MAX_SIZE = 3000
_MINMAPQUALI = 25
_MIN_COVERAGE = 3
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


def call_output(starts, output_dic, counter_idx,
                max_lst_range=_MAX_SIZE, max_size=_MAX_SIZE):
    if starts:
        starts_key_sort = sorted(starts)
        while max(starts_key_sort) - min(starts_key_sort) > max_lst_range:
            old_pos = starts_key_sort.pop(0)
            old_count = starts.pop(old_pos, None)
            output_dic[old_pos] = []
            for current in starts_key_sort:
                length = current - old_pos
                if length >= max_size:
                    break
                    ## do no know if break or contiune
                if old_count >= _MIN_COVERAGE:
                    counter_idx.count += 1
                    try:
                        output_dic[old_pos].append(length)
                    except KeyError:
                        output_dic[old_pos] = [length]


def update(dic, pos):
    try:  # only True if present in dict
        dic[pos] += 1
    except KeyError:
        dic[pos] = 1


def writetofile(output_dic, f_name):
    ''' dfs '''
    with gzip.open('{}.gz'.format(f_name), 'w') as f_output:
        for k in sorted(output_dic.iterkeys()):
            for v_each in output_dic[k]:
                f_output.write('{}\t{}\n'.format(k, v_each))


def main(argv):
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    output_dic = {}
    starts = {}  # initialize the positions and counts
    ends = {}
    last_tid = -1
    counter_idx = Phaso_count()
    for chrom, start, end in read_bed(args):
        starts = {}
        ends = {}

        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            if last_tid != record.tid:
                last_tid = record.tid
                call_output(starts, output_dic, counter_idx, max_lst_range=0)
                call_output(ends, output_dic, counter_idx, max_lst_range=0)
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
                call_output(present_dic, output_dic, counter_idx)
            if counter_idx.count > _OUTLENGTH:
                writetofile(output_dic, args.out)
                sys.exit()

        call_output(ends, output_dic, counter_idx, max_lst_range=0)
        call_output(starts, output_dic, counter_idx, max_lst_range=0)
    writetofile(output_dic, args.out)
    samfile.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
