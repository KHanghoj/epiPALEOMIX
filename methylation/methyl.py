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
import fileinput

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
    headers = 'positions not_changed C_to_T C_to_other sum'.split()
    sites = 0
    f = open('dingdong.txt', 'w')
    f.write('\t'.join(headers)+'\n')

    for pileupcolumn in samfile.pileup(args.chrom, args.start, args.end,
                                       truncate=True):

        sites += 1
        posi = pileupcolumn.pos-_FASTAIDX-1
        if not 'CG' in fastaread[posi:posi+2]:
            continue
        not_changed = 0
        C_to_T = 0
        C_to_other = 0
        tot = 0

        for x in pileupcolumn.pileups:

            if x.alignment.is_reverse:

                if (x.qpos == x.alignment.alen-1) and  \
                        x.alignment.seq[x.qpos-1] == 'C' \
                        and (not x.indel):

                    tot += 1
                    if x.alignment.seq[x.qpos] == 'A':
                        C_to_T += 1
                    if x.alignment.seq[x.qpos] == 'G':
                        not_changed += 1
                    if x.alignment.seq[x.qpos] in 'CT':
                        C_to_other += 1

            else:
                if x.qpos == 1 and  \
                        x.alignment.seq[x.qpos] == 'G' \
                        and (not x.indel):

                    tot += 1

                    if x.alignment.seq[x.qpos-1] == 'T':
                        C_to_T += 1
                    if x.alignment.seq[x.qpos-1] == 'C':
                        not_changed += 1
                    if x.alignment.seq[x.qpos-1] in 'AG':
                        C_to_other += 1

        print(pileupcolumn.pos+1, not_changed,
              C_to_T, C_to_other, tot, file=f, sep='\t')
    print(sites)

    f.close()
    samfile.close()
    fasta.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))










