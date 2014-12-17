#!/opt/local/bin/python
'''  Object: call the nucleosome center if exists and return key=POSITION,
value=DEPTH. It extends the values around the centre of the nucleosome if
identical maximal depth
Example:
test.bam --chrom 22 --start 18500000 --end 18513769 --out dingdong
Thi is with and empty '\n' bed file.
nucleomap_bedformat.py test.bam test.bed --chrom 22 --start 20000000 --end 20200000 --out tra12
'''
from __future__ import print_function

#from fileinput import inpua
import sys
import pysam
import math
import argparse

#### Constants:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)


def call_max(mainwind):
    ''' docstring '''
    call = {}  # this is the return dic with position and depth
    maxdepth = max(mainwind)
    if maxdepth <= 2:  # do not take maxseqdep of 2 or lower into account
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

        # for i in xrange(center_index - 1, -1, -1):
        #     if mainwind[i] == maxdepth):
        #         call[i] = mainwind[i]
        #     else:
        #         break

        # for i in xrange(center_index + 1, len(mainwind)):
        #     if mainwind[i] == maxdepth):
        #         call[i] = mainwind[i]
        #     else:
        #         break
    else:
        return None  # not callable. i.e the centre is not maximal

    return call


def call_flanking(flankwindow):
    ''' docstring '''
    return sum(flankwindow)/len(flankwindow)


def shift_window(pileupcolumn, windows, positions, last_pos, s_depth):
    ''' docstring '''
    delta_pos = pileupcolumn.pos - last_pos
    if delta_pos >= _TOTAL_WIN_LENGTH:
        windows[:] = [0] * _TOTAL_WIN_LENGTH
        positions[:] = [(pileupcolumn.tid, last_pos + pos)
                        for pos in xrange(pileupcolumn.pos - _TOTAL_WIN_LENGTH,
                                          pileupcolumn.pos)]
    else:
        for pos in xrange(1, delta_pos):
            windows.append(0)
            positions.append((pileupcolumn.tid, last_pos + pos))
    windows.append(s_depth)  # the actual windows with pileup
    positions.append((pileupcolumn.tid, pileupcolumn.pos))  # positions

    while len(windows) > _TOTAL_WIN_LENGTH:
        windows.pop(0)
        positions.pop(0)


def call_window(windows, positions, last_result, samfile, outname):
    ''' docstring '''
    calls_center = call_max(windows[_NEIGHBOR + _OFFSET: _NEIGHBOR +
                                    _OFFSET+_SIZE])

    if calls_center:  # checks if the center of nuclesome == maxdepth
        calls_spacerL = call_flanking(windows[:_NEIGHBOR])
        calls_spacerR = call_flanking(windows[-_NEIGHBOR:])

        if not (calls_spacerL is None or calls_spacerR is None):
            mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)

            value = next(calls_center.itervalues())
            score = math.log(float(value) / mean_spacer)
            # math.log by standard is natural
            position_start = positions[_NEIGHBOR + _OFFSET + min(calls_center)]
            position_end = positions[_NEIGHBOR + _OFFSET + max(calls_center)]

            # NOTE: Use getrname to get name of chromosome / contig / etc.
            chrom = samfile.getrname(position_start[0])
            result = [chrom, position_start[1]+1, position_end[1]+1,
                      value, score]

            if result[:3] == last_result[:3]:  # to avoid printing duplicates
                if last_result[4] > score:
                    result[4] = last_result[4]
            else:  # result[:3] != last_result[:3]:
                if not last_result[0] is 0:
                    print(*last_result, file=outname, sep='\t')
            # return result
            del last_result[:]
            last_result.extend(result)


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


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    f_output = open(args.out, 'w')
    windows = []
    last_tid = -1
    last_pos = -1
    sizeofoutputfile = 5
    last_result = [0]*sizeofoutputfile
    positions = []  # the chromosomal position of the nuclesome
    samfile = pysam.Samfile(args.bam, "rb")
    for chrom, start, end in read_bed(args):
        last_tid = -1
        last_pos = -1
        windows = []
        last_result = [0]*sizeofoutputfile
        for pileupcolumn in samfile.pileup(chrom, start, end):
            if pileupcolumn.tid != last_tid:
                last_tid = pileupcolumn.tid
                last_pos = -1
                windows = []
                last_result = [0]*sizeofoutputfile

            s_depth = int(pileupcolumn.n)
            shift_window(pileupcolumn, windows, positions,
                         last_pos, s_depth)
            last_pos = pileupcolumn.pos
            if len(windows) == _TOTAL_WIN_LENGTH:
                call_window(windows, positions,
                            last_result, samfile, f_output)

    print(*last_result, file=f_output, sep='\t')
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
