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


def get_size(dic_start, dic_end, samfile, chrom, f_output):
    start_pos = min(dic_start)
    print(start_pos)
    startstrand = dic_start[min(dic_start)]
    end_pos = max([k for k, v in dic_end.items() if v != startstrand])
    end_strand = dic_end[end_pos]
    length = end_pos-start_pos
    print(samfile.getrname(chrom), start_pos, end_pos,
          length, startstrand, end_strand, file=f_output, sep='\t')


def main(argv):
    ''' docstring '''
    chrom = ''
    currentend = -1
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")

    headers = 'chrom\tstart\tend\tlength\tstarts_reverse'

    f_output = open(args.out, 'w')  # the output file
    f_output.write(headers+'\n')
    dic_start = {}
    dic_end = {}
    for record in samfile.fetch(args.chrom, args.start, args.end):
        increment = record.pos + record.qend
        
        if record.tid != chrom:  # new chromosome
            ## print current result IMPORTANT
            chrom = record.tid
            currentend = increment
            dic_start.clear()
            dic_end.clear()
        print(currentend, record.pos)
        if currentend >= record.pos:  # if new read overlaps former
            currentend = increment
            dic_start[record.pos] = record.is_reverse
            dic_end[currentend] = record.is_reverse
            print(dic_start)
            print(dic_end)
        else:
            get_size(dic_start, dic_end, samfile, chrom, f_output)
            # start_pos = min(dic_start)
            # print(start_pos)
            # startstrand = dic_start[min(dic_start)]
            # end_pos = max([k for k, v in dic_end.items() if v != startstrand])
            # end_strand = dic_end[end_pos]
            # length = end_pos-start_pos
            # print(samfile.getrname(chrom), start_pos, end_pos,
            #       length, startstrand, end_strand, file=f_output, sep='\t')

            dic_start.clear()
            dic_end.clear()
            chrom = ''
    get_size(dic_start, dic_end, samfile, chrom, f_output)
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
