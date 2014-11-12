#!/opt/local/bin/python
'''  hahaha doctstring '''
from __future__ import print_function

#from fileinput import inpua
import sys
import pysam
import math
import collections
import random


# NOTE: Global constants are normally written UPPERCASE, private varibles
# starting with single underscore (_SIZE, _OFFSET, etc.)

#### Constants:
size = 147 # tpenhe window
offset = 12 # offset between spacer and nucleosomal DNA
neighbor = 25  #size of flanking regions to be considered for the log-odd ration score 
_NEIGHBOR = 25
# NOTE: Moved from below
total_win_len = size+(2*offset)+(2*neighbor)

'''  Object: call the nucleosome center if exists and return
key=POSITION, value=DEPTH. It extends the values around the
centre of the nucleosome if identical maximal depth '''


def call_max(sliceT):
    ''' sadfafasdfds'''
    call = {}  # this is the return dic with position and depth
    maxdepth = max(sliceT)
    if maxdepth <= 2:  # do not take maxseqdep of 2 or lower into account
        # NOTE: Return None instead of setting 'NA' = 1
        return None
    center_index = (len(sliceT)-1)/2

    if sliceT[center_index] == maxdepth:
        call[center_index] = sliceT[center_index]
        i = 1
        # NOTE: Use 'and' instead of '&' unless non-lazy evaluation is needed
        while (sliceT[(center_index - i)] == maxdepth) & (i < center_index):
            call[(center_index-i)] = sliceT[(center_index-i)] ## extends 5'
            i+=1

        i=1
        while (sliceT[(center_index + i)] == maxdepth) & (i < center_index):
            call[(center_index+i)] = sliceT[(center_index+i)] ## extends 3'
            i+=1

        # for i in xrange(center_index - 1, -1, -1):
        #     if sliceT[i] == maxdepth):
        #         call[i] = sliceT[i]
        #     else:
        #         break

        # for i in xrange(center_index + 1, len(sliceT)):
        #     if sliceT[i] == maxdepth):
        #         call[i] = sliceT[i]
        #     else:
        #         break
    else:
        return None  # not callable. i.e the centre is not maximal

    return call


def call_flanking(sliceT):
    maxdepth = max(sliceT)
    center_index = (len(sliceT)-1)/2
    if sliceT[center_index] == maxdepth:
        return sliceT[center_index]


def shift_window(pileupcolumn, windows, positions, last_pos, s_depth): 
    delta_pos = pileupcolumn.pos - last_pos
    if delta_pos >= total_win_len:
        windows[:] = [0] * total_win_len
        positions = [(pileupcolumn.tid, last_pos + pos)
                     for pos in xrange(pileupcolumn.pos - total_win_len,
                                       pileupcolumn.pos)]
    else:
        for pos in xrange(1, delta_pos):
            windows.append(0)
            positions.append((pileupcolumn.tid, last_pos + pos))

    windows.append(s_depth) # the actual windows with pileup
    positions.append((pileupcolumn.tid, pileupcolumn.pos)) # positions

    while len(windows) > total_win_len:
        windows.pop(0)
        positions.pop(0)


def call_window(windows, positions, last_result, samfile):
    calls_center = call_max(windows[neighbor+offset:neighbor+
        offset+size])

    if calls_center: # checks if the center of nuclesome == maxdepth
        calls_spacerL = call_flanking(windows[:neighbor])
        calls_spacerR = call_flanking(windows[-neighbor:])

        if not (calls_spacerL is None or calls_spacerR is None):
            mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)

            value = next(calls_center.itervalues())
            score = math.log(float(value) / mean_spacer)
            position_start = positions[neighbor + offset + min(calls_center)]
            position_end = positions[neighbor + offset + max(calls_center)]

            # NOTE: Use getrname to get name of chromosome / contig / etc.
            chrom = samfile.getrname(position_start[0])
            result = (chrom, position_start[1], position_end[1],
                      value, score)
            if result != last_result:
                print(*result)

            return result

### here starts the script:
# NOTE: These variables belong in a function:
sites = 0
windows = []   # the sum of 25+12+147+12+25=221
last_tid = -1
last_pos = -1
last_result = None
positions = []  # the chromosomal position of the nuclesome
basecomposi = []  # single nucleotide base composition
posi = []
cnt = collections.Counter()  # this is for counting bases in basecomposi


# NOTE: Check number of arguments


f = argv[1]
samfile = pysam.Samfile(f, "rb")
# result = list(x.n for x in samfile.pileup())


for pileupcolumn in samfile.pileup():   #'22',16050360,16057100):
    if pileupcolumn.tid != last_tid:
        last_tid = pileupcolumn.tid
        last_pos = -1
        windows = []

    # NOTE: Do not include positions that are indels
    bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
             if not x.indel]
    # NOTE: Don't assume that base 0 is the consensus
    # get the consensus nucleotide at each site
    #  TODO: Handle deletions by checking if bases is empty
    #  base = random.choice(bases)
    s_depth = len(bases)  # get the depth

    shift_window(pileupcolumn, windows, positions, last_pos, s_depth)
    last_pos = pileupcolumn.pos

    if len(windows) == total_win_len:
        last_result = call_window(windows, positions, last_result, samfile)

samfile.close()

import argparse


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('name', help="...")
    parser.add_argument('--param', help="...")
    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)
    print(args.name)
    print(args.param)

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))



# from collections import defaultdict
# def C(s):
#     dinucleot = 2
#     d = defaultdict(int)
#     # d = {}
#     for i in xrange(len(s)-dinucleot-1):
#         d[s[i:i+dinucleot]] += 1
#     return d

# print(C('ACGTAGCTACGACTCGATGCATGTAGCTAGCTAGCTACGTCAGTAGCTGACTGATCGATCGTCGTAGCTGACTGATCGATCGATCGAATTCCGG'))
