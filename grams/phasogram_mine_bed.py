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
from collections import defaultdict

_MIN_DEPTH = 5
_MAX_SIZE = 1000
_MINMAPQUALI = 30
_NEXTNUC = 10


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


def get_distance(start, end, score, output_dic):
    ''' docstring '''
    reads_for_removal = []
    last_element = len(end)-1
    for idx in range(last_element-1):  # do not count itself
        if start[idx] <= end[last_element]:
            score[idx] += 1
            score[last_element] += 1
        if abs(end[idx]-start[last_element]) > _MAX_SIZE:
            reads_for_removal.append(idx)  # get the index

    for idx_remove in reads_for_removal:
        if score[idx_remove] >= _MIN_DEPTH:
            for idx_list in range(last_element):
                if not idx_remove == idx_list:
                    var = abs(start[idx_remove]-start[idx_list])
                    if var >= _NEXTNUC:
                        output_dic[var] += 1
                        # this creates much smaller output files.
                        # f_output.write(str(var)+'\n')

    # removal:
    for idx in sorted(reads_for_removal, reverse=True):
        start.pop(idx)
        end.pop(idx)
        score.pop(idx)


def writetofile(output_dic, f_output):
    for key, value in output_dic.iteritems():
        f_output.write('{}\t{}\n'.format(key, value))


def main(argv):
    ''' docstring '''
    chrom_check = ''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    output_dic = defaultdict(int)
    try:
        remove(args.out)
        f_output = open(args.out, 'a')  # the output file
    except OSError:
        f_output = open(args.out, 'a')  # the output file
    start_plus, start_minus = [], []
    end_plus, end_minus = [], []
    score_plus, score_minus = [], []

    with open(args.bed, 'r') as myfile:
        for line in myfile.readlines():
            input_line = (line.rstrip('\n')).split('\t')[:3]
            chrom = input_line.pop(0).replace('chr', '')
            start = int(input_line.pop(0))
            end = int(input_line.pop(0))
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

                    get_distance(start_minus, end_minus,
                                 score_minus, output_dic)
                else:
                    get_append(start_plus, end_plus, score_plus,
                               record.pos, record.aend)
                    get_distance(start_plus, end_plus, score_plus, output_dic)
    writetofile(output_dic, f_output)
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
