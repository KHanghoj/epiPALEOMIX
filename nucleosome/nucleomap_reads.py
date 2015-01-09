#!/opt/local/bin/python
from __future__ import print_function

#from fileinput import inpua
import sys
import pysam
import math
import argparse
from collections import deque
from itertools import islice
#### Constants:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_MAXLEN = int(1e3)
_MINDEPTH = 2


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
        # NOTE: Use 'and' instead of '&' unless non-lazy evaluation is needed
        while (mainwind[(center_index - i)] == maxdepth) and (i < center_index):
            call[(center_index-i)] = mainwind[(center_index-i)]  # extends 5'
            i += 1
        i = 1
        while (mainwind[(center_index + i)] == maxdepth) and (i < center_index):
            call[(center_index+i)] = mainwind[(center_index+i)]  # extends 3'
            i += 1
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


def call_window(depths_deque, positions_deque, chrom, output_dic):
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
                # math.log by default is natural
                win_positions = list(win_positions)
                position_start = win_positions[_NEIGHBOR +
                                               _OFFSET + min(calls_center)]
                position_end = win_positions[_NEIGHBOR +
                                             _OFFSET + max(calls_center)]

                key = '{}_{}_{}'.format(chrom, position_start,
                                        position_end)
                if key in output_dic.keys():
                    # if new score is great than old score
                    # at same positions then replace.
                    if output_dic[key][1] < score:
                        output_dic[key] = [value, score]
                else:
                    output_dic[key] = [value, score]


def fun(item):
    try:
        return int(item)
    except ValueError:
        return str(item)


def writetofile(output_dic, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    key_values = iter(sorted((map(fun, key.split('_'))+value for key,
                      value in output_dic.iteritems())))
    fmt = '{0}\t{1}\t{2}\t{3}\t{4}\n'
    for dat in key_values:
        f_output.write(fmt.format(*dat))


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    last_tid = ''
    last_pos = -1
    output_dic = {}
    samfile = pysam.Samfile(args.bam, "rb")

    for chrom, start, end in read_bed(args):
        last_tid = ''
        last_pos = -1

        depths_deque = deque(maxlen=_MAXLEN)
        positions_deque = deque(maxlen=_MAXLEN)
        depths_deque = []
        positions_deque = []
        for record in samfile.fetch(chrom, start, end):
            if record.tid != last_tid:
                #########
                ######### Remember to call window before new chrom.
                #########
                last_tid = record.tid
                last_pos = record.pos
                # depths_deque.clear()
                depths_deque.extend([0]*_MAXLEN)
                # positions_deque.clear()
                positions_deque.extend((pos+1) for pos in xrange(record.pos,
                                       record.pos+_MAXLEN))
                deque_idx = 0

            delta_pos = record.pos - last_pos
            deque_idx += delta_pos
            if deque_idx + record.alen >= _MAXLEN:
                chrom = samfile.getrname(record.tid)
                call_window(depths_deque, positions_deque, chrom, output_dic)
                ## run the scoring through the depths_deque
                # depths_deque.clear()
                depths_deque.extend([0]*(_MAXLEN-_TOTAL_WIN_LENGTH))
                # positions_deque.clear()
                positions_deque.extend((pos+1) for pos in xrange(record.pos,
                                       record.pos+_MAXLEN-_TOTAL_WIN_LENGTH))

                while depths_deque < _MAXLEN:
                    depths_deque.pop(0)
                    positions_deque.pop(0)
                # positions_deque.extend((record.tid,
                #                        pos+1) for pos in xrange(record.pos,
                #                        record.pos+_MAXLEN-_TOTAL_WIN_LENGTH))
                deque_idx = _TOTAL_WIN_LENGTH
                # delta_pos = _TOTAL_WIN_LENGTH
            # deque_idx += delta_pos
            # if deque_idx+record.alen >= _MAXLEN:
            #     # if deque is too long, need to more depths_deque
            #     ## run the scoring through the depths_deque
            #     ## extend until all except the last window. _TOTAL_WIN_LENGTH
            #     depths_deque.extend([0]*(_MAXLEN-_TOTAL_WIN_LENGTH))
            #     deque_idx = 0
            update_depth(depths_deque, record, deque_idx)

            last_pos = record.pos
            last_tid = record.tid
            # print(deque_idx)
    # print(list(islice(depths_deque, 16, 20)))
    # print(list(islice(positions_deque, 16, 20)))
    # print(depths_deque)
    # print(positions_deque)
    # print(zip(depths_deque,positions_deque))
    writetofile(output_dic, args.out)
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
