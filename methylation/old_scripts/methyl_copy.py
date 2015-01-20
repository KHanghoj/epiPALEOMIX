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

    def __init__(self, filename, seq_len=1e6):
        self._fasta = pysam.Fastafile(filename)
        self._seq_len = int(seq_len)
        self._last_chrom = None
        self._fasta_str = None
        self._last_start = None
        self._actual_pos = None
        self._end = None

    def fetch_string(self, chrom, start, nbases=2):
        ''' docstring '''
        if self._last_chrom != chrom or (start-self._last_start) >= \
                self._seq_len or start >= self._end - nbases:

            self._end = start + self._seq_len
            self._fasta_str = self._fasta.fetch(chrom,
                                                start=start, end=self._end)
            self._last_start = start
            self._last_chrom = chrom
        self._actual_pos = start-self._last_start
        return self._fasta_str[self._actual_pos:self._actual_pos+nbases]

    def closefile(self):
        ''' docstring '''
        return self._fasta.close()


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
    top = dic_base_forward.get('T', 0)+tempdic_minus.get('A', 0)
    lower = top+dic_base_forward.get('C', 0)+tempdic_minus.get('G', 0)
    ms_value = top/float(lower)
    print(chrom, last_pos+1, ms_value, file=output, sep='\t')


def call_minus_ms(chrom, last_pos, dic_lastpos, output):
    ''' docstring '''
    for keys in sorted(dic_lastpos.keys()):
        if last_pos > keys:  # this is not necessary
            call_ms(chrom, keys, dic_lastpos, {}, output)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = Cache(args.fastafile)
    # MS: Consider using WITH statements
    f_output = open(args.out, 'w')  # the output file
    chrom = None
    dic_lastpos = defaultdict(lambda: defaultdict(int))

    dic_base_forward = defaultdict(int)

    for record in samfile.fetch(args.chrom, args.start, args.end):
        read_sequence = record.seq
        if record.tid != chrom:  # new chromosome or first record
            chrom = record.tid
            last_pos = -1
            present_chrom = samfile.getrname(record.tid)

        # Call minus strand Ms with no plus strand information
        if dic_lastpos and max(dic_lastpos.keys()) < last_pos:
            call_minus_ms(present_chrom,
                          last_pos, dic_lastpos, f_output)

        if record.is_reverse:  # the minus strand
            if read_sequence[-2:] in _MINUS_STRAND_BASES:  # last two bases ok
                pos = record.aend-2  # end of read minus 2 bases
                if 'CG' in fasta.fetch_string(present_chrom, pos):
                    cigar_op, cigar_len = record.cigar[-1]
                    if (cigar_op == 0) and (cigar_len >= 2):
                        dic_lastpos[record.aend-2][read_sequence[-1]] += 1

        else:  # this is for the forward strand
            if read_sequence[:2] in _PLUS_STRAND_BASES:  # first two bases ok
                if 'CG' in fasta.fetch_string(present_chrom, record.pos):
                    cigar_op, cigar_len = record.cigar[0]
                    if cigar_op == 0 and cigar_len >= 2:
                        if record.pos != last_pos and last_pos != -1:
                            call_ms(present_chrom, last_pos,
                                    dic_lastpos, dic_base_forward, f_output)
                            dic_base_forward.clear()
                        dic_base_forward[read_sequence[0]] += 1
                        last_pos = record.pos

    if dic_base_forward:  # this checks if dic contains anything
        call_ms(present_chrom, last_pos,
                dic_lastpos, dic_base_forward, f_output)
    if dic_lastpos:
        call_minus_ms(present_chrom,
                      last_pos, dic_lastpos, f_output)

    f_output.close()
    samfile.close()
    fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))