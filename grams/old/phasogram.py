#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse

_MIN_DEPTH = 5
_MAX_SIZE = 1000
_MINMAPQUALI = 30
_NEXTNUC = 0


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_phasogram.txt')
    return parser.parse_args(argv)


def empty_lists(*arg):
    for x in arg:
        del x[:]


def get_append(start, end, score, record_begin, record_end):
    start.append(record_begin)
    end.append(record_end)
    score.append(0)


def get_distance(start, end, score, f_output):
    reads_for_removal = []
    last_element = len(end)-1
    for idx in range(last_element-1):  # do not count itself
        if start[idx] <= end[last_element]:
            score[idx] += 1
            score[last_element] += 1
        if abs(end[idx]-start[last_element]) > _MAX_SIZE:
            reads_for_removal.append(idx)  # get the index

    for idx in reads_for_removal:
        if score[idx] >= _MIN_DEPTH:
            for idx_list in range(last_element):
                if not idx == idx_list:
                    var = abs(start[idx]-start[idx_list])
                    if var >= _NEXTNUC:
                        print(var, file=f_output)
    ### removal:
    for idx in sorted(reads_for_removal, reverse=True):
        start.pop(idx)
        end.pop(idx)
        score.pop(idx)


def main(argv):
    ''' docstring '''
    chrom = ''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    f_output = open(args.out, 'w')  # the output file
    start_plus, start_minus = [], []
    end_plus, end_minus = [], []
    score_plus, score_minus = [], []

    for record in samfile.fetch(args.chrom, args.start, args.end):
        if record.mapq < _MINMAPQUALI:
            break  # do not analyze low quality reads
        if chrom != record.tid:
            chrom = record.tid
            empty_lists(start_plus, start_minus, end_plus, end_minus,
                        score_plus, score_minus)  # empty all lists)

        if record.is_reverse:  # minus strand
            get_append(start_minus, end_minus, score_minus,
                       record.aend, record.pos)  # begin and end swapped
            get_distance(start_minus, end_minus, score_minus, f_output)
        else:  # plus strand
            get_append(start_plus, end_plus, score_plus,
                       record.pos, record.aend)
            get_distance(start_plus, end_plus, score_plus, f_output)

    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
