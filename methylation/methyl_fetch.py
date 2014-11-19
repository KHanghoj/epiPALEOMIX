#!/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
 python ~/research/projects/epiomix/methylation/methyl_fetch.py chrome.fa test.bam --chrom 22 --start 18100000 --end 20100000 --out new.txt
 Rscript -e "a=read.table('new.txt') ;summary(a)"
        # Profiling af python script:
        # $ python -m cProfile -s cumulative script.py

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
    parser.add_argument('--out', help='...', default='out_mapmethyl_fetch.txt')
    return parser.parse_args(argv)


def get_ms(chrom, last_pos, dic_lastpos, dic_base_forward, output):
    tempdic_minus = dic_lastpos.pop(last_pos, {})
    top = dic_base_forward['T']+tempdic_minus.get('A', 0)
    lower = top+dic_base_forward['C']+tempdic_minus.get('G', 0)
    M_value = top/float(lower)
    dic_base_forward.clear()
    print(chrom, last_pos+1, M_value,
          file=output, sep='\t')


def get_minus_ms(chrom, last_pos, dic_lastpos, output):
    for keys in sorted(dic_lastpos.keys()):
        if last_pos > keys:
            dic_temp = dic_lastpos.pop(keys, {})
            top = dic_temp.get('A', 0)
            lower = top + dic_temp.get('G', 0)
            M_value = top/float(lower)
            print(chrom, keys+1, M_value,
                  file=output, sep='\t')


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = pysam.Fastafile(args.fastafile)
    f_output = open(args.out, 'w')  # the output file
    last_pos = -1
    chrom = ''

    dic_lastpos = defaultdict(lambda: defaultdict(int))
    dic_base_forward = defaultdict(int)

    # mybeds = read.beds(...)
    # handle = pysam.samfile(...)
    # for bed in mybeds:
    #    call_Ms(handle, bed)

# compare perforance with former script
# see if i can avoid the keys thing
# i get    2578 lines in new.txt when len(dic_lastpos.keys()) >= 2
# i get     2513 lines in new.txt len(dic_lastpos.keys()) >= 200


    for record in samfile.fetch(args.chrom, args.start, args.end):
        read_sequence = record.seq
        read_cigar = record.cigar

        if len(dic_lastpos.keys()) > 0 and (max(dic_lastpos.keys()) < last_pos):
            get_minus_ms(samfile.getrname(record.tid),
                         last_pos, dic_lastpos, f_output)

        if record.tid != chrom:  # new chromosome
            chrom = record.tid
            last_pos = -1

        if record.is_reverse:  # the minus strand
            if read_sequence[-2:] in _MINUS_STRAND_BASES:  # last two bases ok
                fetch_posi = record.aend-2
                if 'CG' in fasta.fetch(samfile.getrname(record.tid),
                                       start=fetch_posi, end=fetch_posi+2):
                    cigar_op, cigar_len = read_cigar[-1]
                    if (cigar_op == 0) and (cigar_len >= 2):
                        dic_lastpos[record.aend-2][read_sequence[-1]] += 1

        else:  # this is for the forward strand
            if read_sequence[:2] in _PLUS_STRAND_BASES:  # first two bases ok
                fetch_posi = record.pos
                if 'CG' in fasta.fetch(samfile.getrname(record.tid),
                                       start=fetch_posi, end=fetch_posi+2):

                    cigar_op, cigar_len = read_cigar[0]
                    if (cigar_op) == 0 and (cigar_len) >= 2:

                        if record.pos != last_pos and last_pos != -1:
                            get_ms(samfile.getrname(record.tid), last_pos,
                                   dic_lastpos, dic_base_forward, f_output)
                        dic_base_forward[read_sequence[0]] += 1
                        last_pos = record.pos
    if len(dic_base_forward.keys()) > 0 or len(last_pos.keys()) > 0:
        get_ms(samfile.getrname(record.tid), last_pos,
               dic_lastpos, dic_base_forward, f_output)

    f_output.close()
    samfile.close()
    fasta.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
