#!/opt/local/bin/python
from __future__ import print_function

#from fileinput import inpua
import sys
import pysam
import math
import argparse
from collections import deque
from itertools import islice
from os import remove
#### Constants:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_MAXLEN = int(1e3)
# _MINDEPTH = 2
_MINDEPTH = 4


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", default=None)
    parser.add_argument('--end', help="...", default=None)
    parser.add_argument('--out', help='...', default='out_nucleomap.txt')
    return parser.parse_args(argv)


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = line.rstrip('\n').split('\t')
                chrom = input_line.pop(0).replace('chr', '')
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        yield (args.chrom, int(args.start), int(args.end))


def update_depth(depths_deque, record, index):
    ## we increment the counts list with all the alli
    for (cigar, count) in record.cigar:
        if cigar in (0, 7, 8):
            for idx in xrange(index, index + count):
                depths_deque[idx] += 1
            index += count
        elif cigar in (2, 3, 6):
            index += count


def call_flanking(flankwindow):
    ''' docstring '''
    return sum(flankwindow)/len(flankwindow)


def call_max(mainwind):
    ''' docstring '''
    call = {}  # this is the return dic with position and depth
    maxdepth = max(mainwind)
    if maxdepth <= _MINDEPTH:
        # NOTE: Return None instead of setting 'NA' = 1
        return None
    center_index = (len(mainwind)-1)/2

    if mainwind[center_index] == maxdepth:
        call[center_index] = mainwind[center_index]
        i = 1
        while (mainwind[(center_index - i)] == maxdepth) and (i < center_index):
            call[(center_index-i)] = mainwind[(center_index-i)]  # extends 5'
            i += 1
        # if wide nucleosome then no good
        if i > 1:
            return None
        i = 1
        while (mainwind[(center_index + i)] == maxdepth) and (i < center_index):
            call[(center_index+i)] = mainwind[(center_index+i)]  # extends 3'
            i += 1
        # if wide nucleosome then no good
        if i > 1:
            return None
    else:
        return None  # not callable. i.e the centre is not maximal

    return call


def window_yield(deque_depth, deque_position):
    iterable_depth = (islice(deque_depth, i,
                      i + _TOTAL_WIN_LENGTH) for i in range(0, _MAXLEN, 1))
    iterable_position = (islice(deque_position, i,
                         i + _TOTAL_WIN_LENGTH) for i in range(0, _MAXLEN, 1))
    while True:
        try:
            output_depth = iterable_depth.next()
            output_position = iterable_position.next()
            yield output_depth, output_position
        except StopIteration:
            break


def call_window(depths_deque, positions_deque, output_dic):
    ''' docstring '''
    for win_depth, win_positions in window_yield(depths_deque, positions_deque):
        win_depth = list(win_depth)
        if not len(win_depth) == _TOTAL_WIN_LENGTH:
            continue
        calls_center = call_max(win_depth[_NEIGHBOR + _OFFSET: _NEIGHBOR +
                                          _OFFSET+_SIZE])
        if calls_center:  # checks if the center of nuclesome == maxdepth
            calls_spacerL = call_flanking(win_depth[:_NEIGHBOR])
            calls_spacerR = call_flanking(win_depth[-_NEIGHBOR:])

            if not (calls_spacerL is None or calls_spacerR is None):
                mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)

                value = next(calls_center.itervalues())
                # because all values are the same
                score = math.log(float(value) / mean_spacer)
                if score > 0:
                    # math.log by default is natural
                    win_positions = list(win_positions)
                    # only start positions is needed
                    # as we only look at 0-wide nucleosomes
                    position_start = win_positions[_NEIGHBOR +
                                                   _OFFSET + min(calls_center)]

                    key = position_start
                    try:
                        output_dic[key] += 1
                    except KeyError:
                        output_dic[key] = 1


def extend_deque(rec_pos, depths_deque,
                 positions_deque, start, deque_idx=None):
    start_pos = rec_pos-_TOTAL_WIN_LENGTH - start
    end_pos = start_pos+_MAXLEN
    if deque_idx is None:
        ## reset all deques
        depths_idx = _MAXLEN
    else:
        ## extend the deques
        depths_idx = deque_idx-_TOTAL_WIN_LENGTH

    depths_deque.extend([0]*(depths_idx))
    positions_deque.extend((pos+1) for pos in xrange(start_pos,
                           end_pos))


def writetofile(output_dic, f_output):
    ''' dfs '''
    ## append to file instead
    mininum = min(output_dic.iterkeys)
    maximum = max(output_dic.iterkeys)
    # key_values = iter(sorted(output_dic.iteritems()))
    fmt = '{0}\t{1}\n'
    # for dat in key_values:
        # f_output.write(fmt.format(*dat))

    for key in xrange(mininum, maximum, 1):
        value = output_dic.get(key, 0)
        f_output.write(fmt.format(key, value))


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    last_tid = ''
    last_pos = -1
    output_dic = {}
    samfile = pysam.Samfile(args.bam, "rb")
    deque_idx = 0
    ## i want to write to file every chrom, to speed up:
    try:
        remove(args.out)
    except OSError:
        pass
    f_output = open(args.out, 'a')
    depths_deque = deque(maxlen=_MAXLEN)
    positions_deque = deque(maxlen=_MAXLEN)

    for chrom, start, end in read_bed(args):
        last_tid = ''
        last_pos = -1

        for record in samfile.fetch(chrom, start, end):
            if record.tid != last_tid:
                if output_dic:
                    call_window(depths_deque, positions_deque, output_dic)
                    depths_deque = deque(maxlen=_MAXLEN)
                    positions_deque = deque(maxlen=_MAXLEN)
                last_tid = record.tid
                last_pos = record.pos
                extend_deque(record.pos, depths_deque, positions_deque, start)
                deque_idx = _TOTAL_WIN_LENGTH

            jump = record.pos - last_pos
            deque_idx += jump

            if jump > _MAXLEN:
                if max(depths_deque) > _MINDEPTH:
                    call_window(depths_deque, positions_deque, output_dic)
                extend_deque(record.pos, depths_deque, positions_deque, start)
                deque_idx = _TOTAL_WIN_LENGTH

            elif deque_idx + record.alen >= _MAXLEN:
                if max(depths_deque) > _MINDEPTH:
                    call_window(depths_deque, positions_deque, output_dic)
                extend_deque(record.pos, depths_deque, positions_deque, start,
                             deque_idx)
                deque_idx = _TOTAL_WIN_LENGTH

            update_depth(depths_deque, record, deque_idx)

            last_pos = record.pos
            last_tid = record.tid
    writetofile(output_dic, f_output)
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
