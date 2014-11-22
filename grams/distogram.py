#!/opt/local/bin/python
'''  Object: To calculate the distograms between mapped reads's start positions
aligning in opposing orientation.
~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --start 16056601 --end 16058000
~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --start 16050500 --end 16050720 --out testdist
16050500 16050720

'''

from __future__ import print_function
import sys
import pysam
import argparse


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_distogram.txt')
    return parser.parse_args(argv)


def call_distance(dic_start, samfile, chrom, f_output):
    ''' docstring '''
    for k in sorted(dic_start.keys()):
        if len(set(dic_start.values())) > 1:
            #   only if data available from both strands
            start_pos = min(dic_start)  # start position
            startstrand = dic_start[min(dic_start)]  # start position strand
            end_pos = max([k for k, v in dic_start.items() if v != startstrand])
            # end_strand = dic_start[end_pos]
            length = end_pos-start_pos

            print(samfile.getrname(chrom), start_pos+1, end_pos+1,
                  length, file=f_output, sep='\t')
            dic_start.pop(start_pos, None)
            # STILL NEED TO CALCULATE LAST EXAMPLE IN VALOEUV 2011
            # need to pop only smallest value and recalculate
    dic_start.clear()


def update_dic(dic_start, beginpos, endpos, strand):
    ''' docstring '''
    if strand:  # minus strand
        dic_start[endpos] = strand
    else:
        dic_start[beginpos] = strand

def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    f_output = open(args.out, 'w')  # the output file
    chrom = ''
    last_end_pos = -1

    dic_start = {}
    # dic_start = {}
    for record in samfile.fetch(args.chrom, args.start, args.end):
        strand = record.is_reverse
        if record.tid != chrom:  # new chromosome or first read
            # call_distance(dic_start, dic_end, samfile, chrom, f_output)
            chrom = record.tid
            last_end_pos = record.aend

        if last_end_pos >= record.pos:  # if new read overlaps former
            last_end_pos = record.aend  # assign new end read
            # if strand:  # minus strand
            #     dic_start[record.aend] = strand
            # else:
            #     dic_start[record.pos] = strand
            update_dic(dic_start, record.pos, last_end_pos, strand)
        else:  # calculate the distance
            call_distance(dic_start, samfile, chrom, f_output)
            last_end_pos = record.aend
            # if strand:  # minus strand
            #     dic_start[record.aend] = strand
            # else:
            #     dic_start[record.pos] = strand
            update_dic(dic_start, record.pos, last_end_pos, strand)


    call_distance(dic_start, samfile, chrom, f_output)
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
