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
    fastaread = fasta.read()  # this is dangerous for the big fasta file.
                              # need to make a loop of some sort or index check

    for pileupcolumn in samfile.pileup(args.chrom, args.start, args.end,
                                       truncate=True):
        posi_rev = pileupcolumn.pos-_FASTAIDX-1
        # if 'CG' in fastaread[posi_rev:posi_rev+2]:
        #     for x in pileupcolumn.pileups:
        #         if x.alignment.is_reverse:
        #             # print('reversed: ', hex(x.alignment.flag))
        #             # print(x.qpos, x.alignment.alen)
        #             if x.qpos == x.alignment.alen-1 and  \
        #                     x.alignment.seq[x.qpos-1] == 'C' and \
        #                     x.alignment.seq[x.qpos] == 'G':
                        # print('reverse')
                        # print(True, pileupcolumn.pos+1)
                        # print(x.alignment.seq[x.qpos-1])
                        # print(x.alignment.seq[x.qpos])

        posi_for = pileupcolumn.pos-_FASTAIDX-1
        if 'CG' in fastaread[posi_for:posi_for+2]:
            for x in pileupcolumn.pileups:
                print(x.qpos==1,x.alignment.seq[x.qpos])
                if x.qpos == 1 and  \
                        x.alignment.seq[x.qpos] == 'G':
                    print('forward')
                    # print(pileupcolumn.pos+1)
                    # print(x.alignment.seq[x.qpos])
                    print(x.alignment.seq[x.qpos-1] == 'C', 'normal')
                    print(x.alignment.seq[x.qpos-1] == 'T', 'deaminated')














    # for pileupcolumn in samfile.pileup(args.chrom, 16051199, 16051199+1,
    #                                    truncate=True):

        # bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
        #          if (x.qpos == 0) and (not x.indel)]

        # bases_end = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
        #              if (x.qpos == len(x.alignment.seq)) and (not x.indel)]
        #             # only the first base in a read
        #          # if not x.indel]  # do not allow insertions or deletions 16051199

        # print([x.qpos for x in pileupcolumn.pileups])  # get position of the base at the given aligned read
        # print([len(x.alignment.seq) for x in pileupcolumn.pileups])  # length of each read aligned to given posisition
        # print([x.alignment.seq[x.qpos] for x in pileupcolumn.pileups])  # get aligned base at position qpos in the read
        # print([a.alignment.cigar for a in pileupcolumn.pileups])  # get cigar for the read
        # print([a.alignment.alen for a in pileupcolumn.pileups])  # alen is the number of bases at for the read
        # print([a.alignment.seq[a.alignment.alen-1] for a in pileupcolumn.pileups])  # the the final base of the read
        # print([a.alignment.flag for a in pileupcolumn.pileups])  # get the flag information. 0x0 is forward 0x10 is reversed
        # print([a.alignment.is_reverse for a in pileupcolumn.pileups])  # whether reversed or not
        # print(pileupcolumn.pos+1)  # get the actual position of the pileup
        # if len(set(bases)) > 1 and fastaread[pileupcolumn.pos-_FASTAIDX] == 'C':
        #     print('bases: ')
        #     print(bases, pileupcolumn.pos)
        #     print('fasta base and position')
        #     print(fastaread[pileupcolumn.pos-_FASTAIDX],
        #           pileupcolumn.pos)
        #     print('BREAKBREAKBREAK')
        #     print(bases_end)
        #     print('BREAKBREAKBREAK')
        #     # '\': line continuation character
        # if len(set(bases_end)) > 1 and  \
        #         fastaread[pileupcolumn.pos-_FASTAIDX] == 'C':
        #     print('BREAKBREAKBREAK')
        #     print('bases_end: ')
        #     print(bases_end, pileupcolumn.pos)
        #     print('fasta base and position')
        #     print(fastaread[pileupcolumn.pos-_FASTAIDX],
        #           pileupcolumn.pos)
        #     # print('head,tail,level')
        #     # print(list(x.is_head for x in pileupcolumn.pileups))  # do not what head,tail,level is? found somethink on samtool. see reading list.
        #     # print(list(x.is_tail for x in pileupcolumn.pileups))
        #     # print(list(x.level for x in pileupcolumn.pileups))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))










