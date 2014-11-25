# !/opt/local/bin/python
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


class Cache(object):
    ''' class doc '''

    def __init__(self, filename):
        self.fasta = pysam.Fastafile(filename)
        self.last_chrom = ''
        self.seq_len = 1e6

    def fetch_string(self, chrom, start):
        if self.last_chrom != chrom or (self.idx+1000) >= self.seq_len:
            self.end = start + self.seq_len
            self.fasta_str = self.fasta.fetch(chrom, start=start, end=self.end)
            self.last_start = start
            self.last_chrom = chrom
            self.idx = 0
        self.idx += start-self.last_start
        self.last_start = start
        return(self.fasta_str[self.idx:self.idx+2])

    def closefile(self):
        return self.fasta.close()


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


def call_ms(chrom, last_pos, dic_lastpos, dic_base_forward, output):
    ''' docstring '''
    tempdic_minus = dic_lastpos.pop(last_pos, {})
    top = dic_base_forward['T']+tempdic_minus.get('A', 0)
    lower = top+dic_base_forward['C']+tempdic_minus.get('G', 0)
    ms_value = top/float(lower)
    dic_base_forward.clear()
    print(chrom, last_pos+1, ms_value, file=output, sep='\t')


def call_minus_ms(chrom, last_pos, dic_lastpos, output):
    ''' docstring '''
    for keys in sorted(dic_lastpos.keys()):
        if last_pos > keys:  # this is not necessary
            dic_temp = dic_lastpos.pop(keys, {})
            top = dic_temp.get('A', 0)
            lower = top + dic_temp.get('G', 0)
            ms_value = top/float(lower)
            print(chrom, keys+1, ms_value, file=output, sep='\t')


def fetchfasta(chrom, presentpos, fasta):
    ''' docstring '''
    end = presentpos + _TEMP_FASTA_LENGTH  # .1 million bases
    return fasta.fetch(chrom, start=presentpos, end=end)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = Cache(args.fastafile)
    f_output = open(args.out, 'w')  # the output file
    chrom = ''
    dic_lastpos = defaultdict(lambda: defaultdict(int))

    dic_base_forward = defaultdict(int)

    for record in samfile.fetch(args.chrom, args.start, args.end):
        present_chrom = samfile.getrname(record.tid)
        read_sequence = record.seq
        read_cigar = record.cigar
        if record.tid != chrom:  # new chromosome or first record
            chrom = record.tid
            last_pos = -1
            # fasta.fetch_string(present_chrom, record.pos)

        # Call minus strand Ms with no plus strand information
        if len(dic_lastpos.keys()) > 0 and max(dic_lastpos.keys()) < last_pos:
            call_minus_ms(samfile.getrname(record.tid),
                          last_pos, dic_lastpos, f_output)

        if record.is_reverse:  # the minus strand
            if read_sequence[-2:] in _MINUS_STRAND_BASES:  # last two bases ok
                pos = record.aend-2  # end of read minus 2 bases
                if 'CG' in fasta.fetch_string(present_chrom, pos):
                    cigar_op, cigar_len = read_cigar[-1]
                    if (cigar_op == 0) and (cigar_len >= 2):
                        dic_lastpos[record.aend-2][read_sequence[-1]] += 1

        else:  # this is for the forward strand
            if read_sequence[:2] in _PLUS_STRAND_BASES:  # first two bases ok
                if 'CG' in fasta.fetch_string(present_chrom, record.pos):
                    cigar_op, cigar_len = read_cigar[0]
                    if cigar_op == 0 and cigar_len >= 2:
                        if record.pos != last_pos and last_pos != -1:
                            call_ms(samfile.getrname(record.tid), last_pos,
                                    dic_lastpos, dic_base_forward, f_output)
                        dic_base_forward[read_sequence[0]] += 1
                        last_pos = record.pos

    # check the absolute final reads
    if len(dic_base_forward.keys()) > 0 or len(dic_lastpos.keys()) > 0:
        call_ms(samfile.getrname(chrom), last_pos,
                dic_lastpos, dic_base_forward, f_output)

    f_output.close()
    samfile.close()
    fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
