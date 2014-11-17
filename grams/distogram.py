#!/opt/local/bin/python
'''  Object: To calculate the distograms between mapped reads's start positions
aligning in opposing orientation.
python ~/research/projects/epiomix/grams/distogram.py test.bam --chrom 22 --end 19100000 --out new.txt
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


def main(argv):
    ''' docstring '''
    chrom = ''
    delta = 25  # if pileupcolumn jump is greater than one, start over
    startstrand_reverse = False
    startstrand_forward = False
    startstrand = ''
    start_pos = -1
    last_pos = -1
    chrom = -1
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")

    headers = 'chrom\tstart\tend\tlength\tforwardstart\treversestart'
    
    f_output = open(args.out, 'w')  # the output file
    f_output.write(headers+'\n')

    for pileupcolumn in samfile.pileup(args.chrom, args.start, args.end,
                                       truncate=True):

        if (pileupcolumn.tid != chrom) or (pileupcolumn.pos - last_pos) > delta:

            # reset all parameters, a new chromosome
            if (last_pos - start_pos) > 0:
                print(samfile.getrname(chrom), start_pos+1, last_pos+1,
                      last_pos-start_pos, startstrand_forward,
                      startstrand_reverse, file=f_output)

            chrom = pileupcolumn.tid
            startstrand = ''
            start_pos = -1
            last_pos = -1
            startstrand_reverse = False
            startstrand_forward = False


        chrom = pileupcolumn.tid
        # get start positions
        # print(startstrand_reverse,pileupcolumn.pos)
        for pileupread in pileupcolumn.pileups:
            # print()
            # print(pileupcolumn.pos, 'pos')
            # print(startstrand, pileupcolumn.pos, 'startstrand')
            # print(pileupread.alignment.is_reverse, 'whether reverse')
            # print(last_pos, 'lastpos')
            # print()
            # here i want to get started
            # print(pileupcolumn.pos)
            # ding =pileupcolumn.pos 
            # if ding >= 16056783-23 and ding <= 16056783:
            #     print(pileupread.alignment.is_reverse,
            #           startstrand_reverse,startstrand_forward,
            #           pileupread.alignment.seq[pileupread.qpos],pileupcolumn.pos)
            #     print(start_pos, last_pos)
            if pileupread.alignment.is_reverse and \
                    (startstrand_reverse is False):  # check strand read
                if last_pos == -1:
                    # startstrand = '-'
                    startstrand_reverse = True
                    start_pos = pileupcolumn.pos
                last_pos = pileupcolumn.pos
                # latest_pos.append(pileupcolumn.pos)
                # positions.append(pileupcolumn.pos)
            elif (pileupread.alignment.is_reverse is False) and \
                 (startstrand_forward is False):
                if last_pos == -1:
                    startstrand_forward = True
                    # startstrand = '+'
                    start_pos = pileupcolumn.pos
                last_pos = pileupcolumn.pos
            # elif int(pileupcolumn.n) >= 1 :
            #     if 
            #     # print(pileupcolumn.n)
            #     last_pos = pileupcolumn.pos


    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
