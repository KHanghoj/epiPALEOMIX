#!/opt/local/bin/python
from __future__ import print_function

#from fileinput import inpua
import sys
import pysam
import math
import argparse
from collections import deque
from itertools import islice, izip, tee
#### Constants:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_MINDEPTH = 4
_DYAD = _OFFSET+_NEIGHBOR+73


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="...")
    parser.add_argument('--out', help='...', default='out_nucleomap.txt')
    return parser.parse_args(argv)


def read_bed(args):
    with open(args.bed, 'r') as myfile:
        for line in myfile.readlines():
            input_line = line.rstrip('\n').split('\t')
            chrom = input_line.pop(0).replace('chr', '')
            start = int(input_line.pop(0))
            end = int(input_line.pop(0))
            yield (chrom, start, end)


def update_depth(depths_lst, record, index):
    ## we increment the counts list with all the alli
    for (cigar, count) in record.cigar:
        if cigar in (0, 7, 8):
            for idx in xrange(index, index + count):
                depths_lst[idx] += 1
            index += count
        elif cigar in (2, 3, 6):
            index += count


def call_flanking(flankwindow):
    ''' docstring '''
    return sum(flankwindow)/len(flankwindow)


def call_max(mainwind):
    ''' docstring '''
    maxdepth = max(mainwind)
    if maxdepth <= _MINDEPTH:
        # NOTE: Return None instead of setting 'NA' = 1
        return None
    center_index = (len(mainwind)-1)/2

    if mainwind[center_index] == maxdepth:
        if (mainwind[(center_index - 1)] == maxdepth or
                mainwind[(center_index + 1)] == maxdepth):
            return None
        return maxdepth


def nwise(deq, n=_TOTAL_WIN_LENGTH):
    ''' produces a list of n=len and increment of 1
        islice ensurest the increment of i, izip zips the first value
        of all lists, the second value of all. Tee simply produces n
        copies of deq '''
    return izip(*(islice(g, i, None)
                  for i, g in enumerate(tee(deq, n))))
    # this is the same output:
    # izip(*(islice(deq, i, None) for i in xrange(n)))


def window_yield(deque_depth, deque_position):
    iterable_depth = nwise(deque_depth)
    iterable_position = nwise(deque_position)
    # iterable_depth = (islice(deque_depth, i,
    #                   i + _TOTAL_WIN_LENGTH) for i in range(0, _MAXLEN, 1))
    # iterable_position = (islice(deque_position, i,
    #                      i + _TOTAL_WIN_LENGTH) for i in range(0, _MAXLEN, 1))
    while True:
        try:
            output_depth = iterable_depth.next()
            output_position = iterable_position.next()
            yield (output_depth, output_position)
        except StopIteration:
            break


def call_window(depths_lst, positions_lst, output_dic):
    ''' docstring '''
    for win_depth, win_positions in window_yield(depths_lst, positions_lst):
        calls_center = call_max(win_depth[_NEIGHBOR + _OFFSET: _NEIGHBOR +
                                          _OFFSET+_SIZE])
        if calls_center:  # checks if the center of nuclesome == maxdepth
            calls_spacerL = call_flanking(win_depth[:_NEIGHBOR])
            calls_spacerR = call_flanking(win_depth[-_NEIGHBOR:])

            if not (calls_spacerL is None or calls_spacerR is None):
                mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)
                # because all values are the same
                score = math.log(float(calls_center) / mean_spacer)
                if score > 0:
                    key = win_positions[_DYAD]
                    try:
                        output_dic[key] += 1
                    except KeyError:
                        output_dic[key] = 1


def writetofile(output_dic, out):
    ''' dfs '''
    with open(out, 'w') as f_output:
        mininum = min(output_dic.iterkeys())
        maximum = max(output_dic.iterkeys())
        fmt = '{0}\t{1}\n'
        for key in xrange(mininum, maximum, 1):
            value = output_dic.get(key, 0)
            f_output.write(fmt.format(key, value))


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    output_dic = {}
    samfile = pysam.Samfile(args.bam, "rb")
    deque_idx = 0
    ## i want to write to file every chrom, to speed up:
    for chrom, start, end in read_bed(args):
        length = end-start
        depths_lst = list([0]*(length+(2*_TOTAL_WIN_LENGTH)))
        positions_lst = xrange(-_TOTAL_WIN_LENGTH, length+_TOTAL_WIN_LENGTH)
        for record in samfile.fetch(chrom, start, end):
            deque_idx = record.pos-start+_TOTAL_WIN_LENGTH
            update_depth(depths_lst, record, deque_idx)
        call_window(depths_lst, positions_lst, output_dic)
    writetofile(output_dic, args.out)
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
