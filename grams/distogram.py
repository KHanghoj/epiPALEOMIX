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
    delta = 1  # if pileupcolumn jump is greater than one, start over
    # reverse = False
    # forward = False
    startstrand = ''
    start_pos = -1
    last_pos = -1
    chrom = -1
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")

    headers = 'chrom start end length '.split()
    f_output = open(args.out, 'w')  # the output file
    f_output.write('\t'.join(headers)+'\n')

    for pileupcolumn in samfile.pileup(args.chrom, args.start, args.end,
                                       truncate=True):

        if (pileupcolumn.tid != chrom) or (pileupcolumn.pos - last_pos) > delta:
            chrom = pileupcolumn.tid
            # reset all parameters, a new chromosome
            # del latest_pos[:]
            if (last_pos - start_pos) != 0:
                print(samfile.getrname(chrom), start_pos+1, last_pos+1,
                      last_pos-start_pos, file=f_output)

            startstrand = ''
            start_pos = -1
            last_pos = -1

        # if (pileupcolumn.pos - last_pos) > delta:
        #     #continue  # to stop counting when out of nucleosome
        #     # here we need to reset paramters as well and print the last known
        #     if (last_pos - start_pos) != 0:
        #         print(samfile.getrname(chrom), start_pos+1, last_pos+1,
        #               last_pos-start_pos, startstrand, file=f_output)
        #     # del latest_pos[:]
        #     startstrand = ''
        #     start_pos = -1
        #     last_pos = -1

            # sys.exit('out of present nucleosome')
        chrom = pileupcolumn.tid
        # get start positions
        for pileupread in pileupcolumn.pileups:
            # print()
            # print(pileupcolumn.pos, 'pos')
            print(startstrand,pileupcolumn.pos, 'startstrand')
            print(pileupread.alignment.is_reverse, 'whether reverse')
            # print(last_pos, 'lastpos')
            # print()
            # here i want to get started
            if pileupread.alignment.is_reverse and \
                    (startstrand is not '-'):  # check strand read
                if last_pos == -1:
                    startstrand = '-'
                    start_pos = pileupcolumn.pos
                last_pos = pileupcolumn.pos
                # latest_pos.append(pileupcolumn.pos)
                # positions.append(pileupcolumn.pos)
            elif (startstrand is not '+') and (pileupread.alignment.is_reverse is False):
                if last_pos == -1:
                    startstrand = '+'
                    start_pos = pileupcolumn.pos
                last_pos = pileupcolumn.pos
            # else:
                # startstrand = ''
                # start_pos = -1
                # last_pos = -1
                # continue
                # latest_pos.append(pileupcolumn.pos)
                # positions.append(pileupcolumn.pos)
            # else: sys.exit('something is wrong')
        # last_pos = pileupcolumn.pos
        # print(pileupcolumn.pos)


        # print(pileupcolumn.pos+1, not_changed,
              # C_to_T, C_to_other, tot_CG, file=f_output, sep='\t')

    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
