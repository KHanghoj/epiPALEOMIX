#!/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
start=19000000, end=20000000 of Fastafile. it just a snippet for testing
python ~/research/projects/epiomix/methylation/methyl.py fasta_snip.fa test.bam --chrom 22 --end 19100000
'''

from __future__ import print_function
import sys
import pysam
import math
import argparse

#### Constants:
_FASTAIDX = 19000000

def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...", default=None)
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_mapmethyl.txt')
    return parser.parse_args(argv)


def getfastafile(file_path):
    # return pysam.Fastafile(file_path)
    return open(file_path, 'r')


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    print(args.fastafile)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = getfastafile(args.fastafile)
    fastaread = fasta.read()

    for pileupcolumn in samfile.pileup(args.chrom, _FASTAIDX, args.end,
                                       truncate=True):

        bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
                 if x.qpos == 0]
                 # if not x.indel]

        if len(set(bases)) > 1 and fastaread[pileupcolumn.pos-_FASTAIDX] == 'C':
            print(bases, pileupcolumn.pos)
            print(fastaread[pileupcolumn.pos-_FASTAIDX],
                  pileupcolumn.pos-_FASTAIDX+_FASTAIDX)
            print()
            print(list(x.is_head for x in pileupcolumn.pileups))
            print(list(x.is_tail for x in pileupcolumn.pileups))
            print(list(x.level for x in pileupcolumn.pileups))
            print(bases)
        # sys.exit('nomore')
    # for base in fastaread:
        # if base == 'A':
            # print(base)
        # else: print('hmm')
    # print(set(fastaread))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
