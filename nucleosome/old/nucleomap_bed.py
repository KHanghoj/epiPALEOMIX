#!/opt/local/bin/python
'''  Object: call the nucleosome center if exists and return key=POSITION,
value=DEPTH. It extends the values around the centre of the nucleosome if
identical maximal depth
Example:
test.bam --chrom 22 --start 18500000 --end 18513769 --out dingdong
Thi is with and empty '\n' bed file.
nucleomap_bedformat.py test.bam test.bed --chrom 22 --start 20000000 --end 20200000 --out tra12
## this is to take i.e. 2e7 numbers:
# start = int(float(input_line.pop(0)))
# end = int(float(input_line.pop(0)))
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
_MINDEPTH = 4
_DYAD = 73


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


def call_flanking(flankwindow):
    ''' docstring '''
    return sum(flankwindow)/len(flankwindow)


def shift_window(pileupcolumn, windows, positions, last_pos, s_depth, start):
    ''' docstring '''
    delta_pos = pileupcolumn.pos - last_pos
    if delta_pos >= _TOTAL_WIN_LENGTH:
        windows[:] = [0] * _TOTAL_WIN_LENGTH
        # positions[:] = [(pileupcolumn.tid, last_pos + pos)
        #                 for pos in xrange(pileupcolumn.pos - _TOTAL_WIN_LENGTH,
        #                                   pileupcolumn.pos)]
        positions[:] = [(pos-start)
                        for pos in xrange(pileupcolumn.pos - _TOTAL_WIN_LENGTH,
                                          pileupcolumn.pos)]
    else:
        for pos in xrange(1, delta_pos):
            windows.append(0)
            positions.append(last_pos + pos-start)

    windows.append(s_depth)  # the actual windows with pileup
    positions.append(pileupcolumn.pos-start)  # positions

    while len(windows) > _TOTAL_WIN_LENGTH:
        windows.pop(0)
        positions.pop(0)


def call_window(windows, positions, samfile, output_dic):
    ''' docstring '''
    calls_center = call_max(windows[_NEIGHBOR + _OFFSET: _NEIGHBOR +
                                    _OFFSET+_SIZE])

    if calls_center:  # checks if the center of nuclesome == maxdepth
        calls_spacerL = call_flanking(windows[:_NEIGHBOR])
        calls_spacerR = call_flanking(windows[-_NEIGHBOR:])

        if not (calls_spacerL is None or calls_spacerR is None):
            mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)
            score = math.log(float(calls_center) / mean_spacer)
            # math.log by default is natural
            if score > 0:
                key = positions[73]
                try:
                    output_dic[key] += 1
                except KeyError:
                    output_dic[key] = 1


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", default=None)
    parser.add_argument('--end', help="...", default=None)
    parser.add_argument('--out', help='...',
                        default='out_nucleomap_bedinput.txt')
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


def writetofile(output_dic, f_output):
    ''' dfs '''
    with open(f_output, 'w') as f:
        mininum = min(output_dic.iterkeys())
        maximum = max(output_dic.iterkeys())
        fmt = '{0}\t{1}\n'
        for key in xrange(mininum, maximum, 1):
            value = output_dic.get(key, 0)
            f.write(fmt.format(key, value))


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    windows = []
    last_tid = ''
    last_pos = -1
    output_dic = {}
    positions = []  # the chromosomal position of the nuclesome
    samfile = pysam.Samfile(args.bam, "rb")
    for chrom, start, end in read_bed(args):
        last_tid = ''
        last_pos = -1
        windows = []
        for pileupcolumn in samfile.pileup(chrom, start, end):
            if pileupcolumn.tid != last_tid:
                # sys.stdout.write(str(pileupcolumn.tid) + '\n')
                last_tid = pileupcolumn.tid
                last_pos = -1
                windows = []

            s_depth = int(pileupcolumn.n)
            shift_window(pileupcolumn, windows, positions,
                         last_pos, s_depth, start)
            last_pos = pileupcolumn.pos
            if len(windows) == _TOTAL_WIN_LENGTH:
                call_window(windows, positions,
                            samfile, output_dic)
    writetofile(output_dic, args.out)
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
