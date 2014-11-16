#!/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
python ~/research/projects/epiomix/methylation/methyl.py fasta_snip.fa test.bam --chrom 22 --end 19100000 --out new.txt
'''

from __future__ import print_function
import sys
import pysam
import argparse


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_mapmethyl.txt')
    return parser.parse_args(argv)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = pysam.Fastafile(args.fastafile)

    headers = 'positions not_changed C_to_T C_to_other sum'.split()
    f_output = open(args.out, 'w')  # the output file
    f_output.write('\t'.join(headers)+'\n')

    for pileupcolumn in samfile.pileup(args.chrom, args.start, args.end,
                                       truncate=True):

        posi_fetch = pileupcolumn.pos-1

        if not 'CG' in fasta.fetch(args.chrom, start=posi_fetch,
                                   end=posi_fetch+2):
                # fetch the two bases of interest from the fasta file
            continue
        not_changed = 0
        C_to_T = 0
        C_to_other = 0
        tot_CG = 0

        for pileupread in pileupcolumn.pileups:

            if pileupread.alignment.is_reverse:  # check strand read

                if (pileupread.qpos == pileupread.alignment.alen-1) and  \
                        pileupread.alignment.seq[pileupread.qpos-1] == 'C' \
                        and (not pileupread.indel):  # read on the rev strand

                    tot_CG += 1
                    if pileupread.alignment.seq[pileupread.qpos] == 'A':
                        C_to_T += 1
                    if pileupread.alignment.seq[pileupread.qpos] == 'G':
                        not_changed += 1
                    if pileupread.alignment.seq[pileupread.qpos] in 'CT':
                        C_to_other += 1

            else:  # this is for the forward strand
                if pileupread.qpos == 1 and  \
                        pileupread.alignment.seq[pileupread.qpos] == 'G' \
                        and (not pileupread.indel):

                    tot_CG += 1

                    if pileupread.alignment.seq[pileupread.qpos-1] == 'T':
                        C_to_T += 1
                    if pileupread.alignment.seq[pileupread.qpos-1] == 'C':
                        not_changed += 1
                    if pileupread.alignment.seq[pileupread.qpos-1] in 'AG':
                        C_to_other += 1

        print(pileupcolumn.pos+1, not_changed,
              C_to_T, C_to_other, tot_CG, file=f_output, sep='\t')

    f_output.close()
    samfile.close()
    fasta.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
