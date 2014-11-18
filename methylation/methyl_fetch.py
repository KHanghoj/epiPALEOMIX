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

    headers = ['positions', 'not_changed', 'C_to_T', 'C_to_other', 'sum']
    f_output = open(args.out, 'w')  # the output file
    f_output.write('\t'.join(headers)+'\n')


    last_pos = -1
    scores = {}
    for record in samfile.fetch(args.chrom, args.start, args.end):
        for gammel_pos in xrange(last_pos, record.pos):
            statistikker = scores.pop(gammel_pos, None)
            if statistikker:
                # print Ms







    for record in samfile.fetch(args.chrom, args.start, args.end):
        # Sekvens = record.seq
        # Alignment = sekvens + cigar


        # Starter med:
          # slut_positioner = {}

        # Hvis + streng

            # 1. Er de første 2 baser TpG eller CpG i sekvensen (.seq)?

            # 2. Er de to baser i FASTA filen CpG i alignment?

            # 3. Er de første 2 baser TpG eller CpG i alignment?

            # 3. Samme position som sidste record, eller ny position?

            # 3a. Samme pos = tæl det som N1 eller N2

            # 3b. Anden pos = Udregne Ms for sidste pos; nulstil tælling; 3a.
            #                 Hent M1/M2 fra slut_positioner:
            #                   slut_positioner.pop(last_pos, [0, 0])

        # Hvis - streng

            # 1. Er de første 2 baser CpG eller CpA i sekvensen (.seq)?

            # 2. Udregne slut-position via cigar-streng

            # 3. Er de sidste to baser i FASTA filen CpG?

            # 4. Er de sidste 2 baser CpA eller CpG?

            # 5. Henter vores [M1, M2] fra slut_positioner med pos som nøgle

            # 6. Incrementer M1 eller M2

        # Når ikke flere records; udregn Ms for sidste position






        # Profiling af python script:
        # $ python -m cProfile -s cumulative script.py













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
                    elif pileupread.alignment.seq[pileupread.qpos] == 'G':
                        not_changed += 1
                    else:
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
