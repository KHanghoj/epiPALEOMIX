#!/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
python ~/research/projects/epiomix/methylation/methyl.py fasta_snip.fa test.bam --chrom 22 --end 19100000 --out new.txt
'''

from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict

_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']


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
    # for record in samfile.fetch(args.chrom, args.start, args.end):
    #     for gammel_pos in xrange(last_pos, record.pos):
    #         statistikker = scores.pop(gammel_pos, None)
    #         if statistikker:
    #             # print Ms


    dic_base_forward = defaultdict(int)
    dic_base_end = defaultdict(int)
    dic_lastpos = defaultdict(int)

    for record in samfile.fetch(args.chrom, args.start, args.end):
        read_sequence = record.seq
        read_cigar = record.cigar
        if record.pos in dic_lastpos:
            dic_lastpos[record.pos]
            pass
        # read_alignment_end = record.aend
        # if len(record.cigar)>1 :
        #     print(read_cigar)
        #     print(read_sequence)
        #     print(len(read_sequence))
        #     print(read_postions)
        #     print(read_alignment_end)

        #     cigar_oper_start = read_cigar[0][0]
        #     cigarlength_start = read_cigar[0][1]


        #     # this i need for minus strand: combine it with the record.aend
        #     cigar_oper_end = read_cigar[-1][0]
        #     cigarlength_end = read_cigar[-1][1]

        #     if(cigar_oper_end is 0 and cigarlength_end >= 2):
        #         print(read_cigar)

        #     # sys.exit(1)
        # if record.pos + len(read_sequence) != record.aend:
        #     print(record.cigar)
        #     print(record.pos)
        #     print(record.alen)
        #     print(record.aend, 'aend')
        #     print(len(read_sequence))
        #     sys.exit()
        # continue
        if record.is_reverse:
            if read_sequence[-2:] in _MINUS_STRAND_BASES:  # last two bases ok
                fetch_posi = record.aend-1
                if 'CG' in fasta.fetch(samfile.getrname(record.tid),
                                       start=fetch_posi, end=fetch_posi+2):
                    if read_cigar[-1][0] == 0 and read_cigar[-1][1] >= 2:
                        dic_lastpos[record.aend][read_sequence[-1]] += 1

                        # if record.aend == last_pos:

                        #     dic_base_end[read_sequence[-1]] += 1
                        # else:
                        #     dic_lastpos[record.aend] = dic_base_end
            # this is the minus strand

        else:
            if read_sequence[0:2] in _PLUS_STRAND_BASES:
                fetch_posi = record.pos+1
                if 'CG' in fasta.fetch(samfile.getrname(record.tid),
                                       start=fetch_posi, end=fetch_posi+2):
                    if read_cigar[0][0] == 0 and read_cigar[0][1] >= 2:
                    # if read_sequence[0:2] in _PLUS_STRAND_BASES and \
                    #         (read_cigar[0][0] == 0 and read_cigar[0][1] > 2):
                        if record.pos == last_pos:
                            dic_base_forward[read_sequence[0]] += 1
                            # count with defaultdic if T (N2) or C (N1)
                            pass
                        else:
                            tempdic_minus = dic_lastpos[last_pos]  # M1 M2
                            top = sum(dic_base_forward['T'], tempdic_minus['A'])
                            lower = sum(top, dic_base_forward['C'], tempdic_minus['G'])
                            M_value = top/float(lower)
                            print(M_value)
                            dic_lastpos.pop(last_pos, None)
                                # returns none of last_pos key is not present
                            dic_base_forward.clear()

        # Sekvens = record.seq
        # Alignment = sekvens + cigar


        last_pos = record.pos
        # Starter med:
          # slut_positioner = {}

        # Hvis + streng

            # 1. Er de foerste 2 baser TpG eller CpG i sekvensen (.seq)?

            # 2. Er de to baser i FASTA filen CpG i alignment?

            # 3. Er de foerste 2 baser TpG eller CpG i alignment?

            # 3. Samme position som sidste record, eller ny position?

            # 3a. Samme pos = tael det som N1 eller N2

            # 3b. Anden pos = Udregne Ms for sidste pos; nulstil taelling; 3a.
            #                 Hent M1/M2 fra slut_positioner:
            #                   slut_positioner.pop(last_pos, [0, 0])

        # Hvis - streng

            # 1. Er de foerste 2 baser CpG eller CpA i sekvensen (.seq)?

            # 2. Udregne slut-position via cigar-streng

            # 3. Er de sidste to baser i FASTA filen CpG?

            # 4. Er de sidste 2 baser CpA eller CpG?

            # 5. Henter vores [M1, M2] fra slut_positioner med pos som noegle

            # 6. Incrementer M1 eller M2

        # Naar ikke flere records; udregn Ms for sidste position






        # Profiling af python script:
        # $ python -m cProfile -s cumulative script.py













        # posi_fetch = pileupcolumn.pos-1

        # if not 'CG' in fasta.fetch(args.chrom, start=posi_fetch,
        #                            end=posi_fetch+2):
        #         # fetch the two bases of interest from the fasta file
        #     continue
        # not_changed = 0
        # C_to_T = 0
        # C_to_other = 0
        # tot_CG = 0

        # for pileupread in pileupcolumn.pileups:
        #     if pileupread.alignment.is_reverse:  # check strand read
        #         if (pileupread.qpos == pileupread.alignment.alen-1) and  \
        #                 pileupread.alignment.seq[pileupread.qpos-1] == 'C' \
        #                 and (not pileupread.indel):  # read on the rev strand

        #             tot_CG += 1
        #             if pileupread.alignment.seq[pileupread.qpos] == 'A':
        #                 C_to_T += 1
        #             elif pileupread.alignment.seq[pileupread.qpos] == 'G':
        #                 not_changed += 1
        #             else:
        #                 C_to_other += 1

        #     else:  # this is for the forward strand
        #         if pileupread.qpos == 1 and  \
        #                 pileupread.alignment.seq[pileupread.qpos] == 'G' \
        #                 and (not pileupread.indel):

        #             tot_CG += 1

        #             if pileupread.alignment.seq[pileupread.qpos-1] == 'T':
        #                 C_to_T += 1
        #             if pileupread.alignment.seq[pileupread.qpos-1] == 'C':
        #                 not_changed += 1
        #             if pileupread.alignment.seq[pileupread.qpos-1] in 'AG':
        #                 C_to_other += 1

        # print(pileupcolumn.pos+1, not_changed,
        #       C_to_T, C_to_other, tot_CG, file=f_output, sep='\t')

    f_output.close()
    samfile.close()
    fasta.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
