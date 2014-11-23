#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
from os import remove

_MIN_DEPTH = 5
_MAX_SIZE = 1000
_MINMAPQUALI = 30
_NEXTNUC = 0


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


def empty_lists(*arg):
    ''' docstring '''
    for a_list in arg:
        del a_list[:]


def get_append(start, end, score, record_begin, record_end):
    ''' docstring '''
    start.append(record_begin)
    end.append(record_end)
    score.append(0)


def get_distance(start, end, score, f_output):
    ''' docstring '''
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
                        f_output.write(str(var)+'\n')
#                       outfile.write("\n".join(itemlist))
                        # this is for lists

    # removal:
    for idx in sorted(reads_for_removal, reverse=True):
        start.pop(idx)
        end.pop(idx)
        score.pop(idx)


def main(argv):
    ''' docstring '''
    chrom_check = ''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    try:
        remove(args.out)
    except OSError:
        pass
    f_output = open(args.out, 'a')  # the output file
    start_plus, start_minus = [], []
    end_plus, end_minus = [], []
    score_plus, score_minus = [], []

    with open(args.bed, 'r') as myfile:
        for line in myfile.readlines():
            input_line = (line.rstrip('\n')).split('\t')[:3]
            #chrom = input_line[0].replace('chr', '')
            chrom = input_line.pop(0).replace('chr', '')
            start = int(input_line.pop(0))
            end = int(input_line.pop(0))
            # print(chrom)
            # print(input_line)
            # start, end = [int(x) for x in input_line[1:]]
            for record in samfile.fetch(chrom, start, end):
                if record.mapq < _MINMAPQUALI:
                    break  # do not analyze low quality reads
                if chrom_check != record.tid:
                    chrom_check = record.tid
                    empty_lists(start_plus, start_minus, end_plus, end_minus,
                                score_plus, score_minus)  # empty all lists

                if record.is_reverse:
                    get_append(start_minus, end_minus, score_minus,
                               record.aend, record.pos)  # begin & end swapped
                    get_distance(start_minus, end_minus, score_minus, f_output)
                else:
                    get_append(start_plus, end_plus, score_plus,
                               record.pos, record.aend)
                    get_distance(start_plus, end_plus, score_plus, f_output)

    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
