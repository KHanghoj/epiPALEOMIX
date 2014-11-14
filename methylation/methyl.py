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
    headers = 'positions not_changed C_to_T C_to_other'.split()
    sites = 0
    f = open('dingdong.txt', 'w')



    for pileupcolumn in samfile.pileup(args.chrom, args.start, args.end,
                                       truncate=True):
        sites+=1
        posi = pileupcolumn.pos-_FASTAIDX-1
        if not 'CG' in fastaread[posi:posi+2]:
            continue
        not_changed = 0
        C_to_T = 0
        C_to_other = 0
        tot =0



        # bases = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups
        #          if ((x.qpos == 0) or(x.qpos == x.alignment.alen-1)) and (not x.indel)]
        # if len(bases)>0: print(bases,pileupcolumn.pos+1)
        for x in pileupcolumn.pileups:
            # note: i want to get alle the positions not deaminated:


            if x.alignment.is_reverse:
                # print('reversed: ', hex(x.alignment.flag))
                # print(x.qpos, x.alignment.alen)
                # if x.qpos == x.alignment.alen-1 and  \
                #         x.alignment.seq[x.qpos-1] == 'C' and \
                #         x.alignment.seq[x.qpos] == 'G':

                if (x.qpos == x.alignment.alen-1) and  \
                        x.alignment.seq[x.qpos-1] == 'C':
                    
                    tot+=1
                    # print('reverse')
                    # print(pileupcolumn.pos+1)
                    # print(x.alignment.seq[x.qpos-1])
                    # print(x.alignment.seq[x.qpos])
                    # print(fastaread[posi:posi+2])
                    if x.alignment.seq[x.qpos] == 'A':  # print('deaminated', pileupcolumn.pos+1)
                        C_to_T += 1
                    if x.alignment.seq[x.qpos] == 'G': not_changed += 1
                    if x.alignment.seq[x.qpos] in 'CT' : C_to_other += 1

            else:
                if x.qpos == 1 and  \
                        x.alignment.seq[x.qpos] == 'G':
                    tot += 1
                    # print('forward')
                    # print(pileupcolumn.pos+1)
                    # print(x.alignment.seq[x.qpos])
                    # print(x.alignment.seq[x.qpos-1] == 'C', 'normal')
                    # print(fastaread[posi:posi+2])
                    if x.alignment.seq[x.qpos-1] == 'T':  # print('deaminated', pileupcolumn.pos+1)
                        C_to_T += 1
                    if x.alignment.seq[x.qpos-1] == 'C': not_changed += 1
                    if x.alignment.seq[x.qpos-1] in 'AG' : C_to_other += 1

        # if not_changed > 0: print('CpG still present: ', not_changed, pileupcolumn.pos+1)
        # if C_to_T > 0: print('Deaminated: ', C_to_T, pileupcolumn.pos+1)
        print(pileupcolumn.pos+1, not_changed,
              C_to_T, C_to_other,tot, file=f, sep='\t')
    print(sites)

    f.close()
    samfile.close()
    fasta.close()










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










