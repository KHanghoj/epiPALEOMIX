# !/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
 python ~/research/projects/epiomix/methylation/methyl_fetch.py chrome.fa test.bam --chrom 22 --start 18100000 --end 20100000 --out new.txt
 Rscript -e "a=read.table('new.txt') ;summary(a)"
 Rscript -e "a=read.table('new.txt');b=sum(a[,1])/sum(a[,2]);print(b)"
        # Profiling af python script:
        # $ python -m cProfile -s cumulative script.py
'''

from __future__ import print_function
import sys
import pysam
import argparse
import re
import gzip
from collections import defaultdict, namedtuple
from itertools import chain

_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']
_BASES_CHECK = 6
_BASES_CHECK = 2
_SKIPBASES = 0  # include the first base of read
_SIZE = _BASES_CHECK-_SKIPBASES-1
_BASE_COMP_ANALYZED = 10


def _build_rev_compl_table():
    ''' this is copied from mikkels simulate.py script
        Remember to ask him '''
    table = ["N"] * 256
    for nt_a, nt_b in zip("AC", "TG"):
        table[ord(nt_a)] = nt_b
        table[ord(nt_b)] = nt_a
    return "".join(table)
_COMPL_TABLE = _build_rev_compl_table()


def _reverse_complement(sequence):
    ''' this is copied from mikkels simulate.py script
    _reverse_complement(insert)'''
    return sequence.translate(_COMPL_TABLE)[::-1]


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

    def fetch_string(self, chrom, start, nbases):
        ''' docstring '''
        if self._last_chrom != chrom or (start-self._last_start) >= \
                self._seq_len or start >= self._end - nbases or \
                start < self._last_start:
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


class Methyl_Level(object):
    """docstring for Methyl_Level"""
    def __init__(self, arg):
        self.arg = arg
        self.fasta = Cache(self.arg.fastafile)
        self.lst_dic_lastpos = self._create_lst_dic()
        self.dic_top = [defaultdict(int) for _ in xrange(_BASES_CHECK-1)]
        self.dic_lower = [defaultdict(int) for _ in xrange(_BASES_CHECK-1)]
        self.indexoflist = xrange(len(self.lst_dic_lastpos))
        self.record, self.start, self.chrom = None, None, None
        self.pat = re.compile('CG')
        self.bed_top, self.bed_low = 0, 0
        self.bed_list = list()
        self.bed_list_app = self.bed_list.append
        self.bed_nametuple = namedtuple('row',
                                        ('chrom start end top lower ratio'))
        self.basecomp_lst_dic = [defaultdict(int) for _
                                 in xrange(_BASE_COMP_ANALYZED)]
        self.basecomp_lst = list()
        self.basecomp_lst_app = self.basecomp_lst.append

    def reset_dict(self, chrom, start, end):
        self.lst_dic_lastpos = self._create_lst_dic()
        self.chrom, self.start, self.end = chrom, start, end
        self.bed_top, self.bed_low = 0, 0

    def _create_lst_dic(self, size=_BASES_CHECK-1):
        return [defaultdict(lambda: defaultdict(int)) for _ in range(size)]

    def _getindexes(self, bases_str):
        ''' returns the 0-based indeces of fasta read'''
        return (m.start() for m in self.pat.finditer(bases_str))

    def update(self, record):
        self.record = record
        if self.record.is_reverse:
            self._reverse_strand()
        else:
            self._forward_strand()

    def call_ms(self):
        for idx in self.indexoflist:
            curr_lst_idx = self.lst_dic_lastpos[idx]
            curr_top = self.dic_top[idx]
            curr_lower = self.dic_lower[idx]
            for key in curr_lst_idx.iterkeys():
                self.keypos = key-self.start
                curr_dic = curr_lst_idx.get(key, {})

                top = (curr_dic.get('T', 0) +
                       curr_dic.get('A', 0))
                lower = (top+curr_dic.get('C', 0) +
                         curr_dic.get('G', 0))

                curr_top[self.keypos] += top
                curr_lower[self.keypos] += lower
                self.bed_top += top
                self.bed_low += lower
        if self.bed_top:
            self.bed_list_app(self.bed_nametuple(self.chrom,
                                                 self.start, self.end,
                                                 self.bed_top, self.bed_low,
                                                 (float(self.bed_top) /
                                                  self.bed_low)))

    def writetofile(self):
        ''' this is done to make sure we have all keys
        present in every read position
        even is key value is zero advantage of defaultdict '''

        with open(self.arg.out, 'w') as f_output:
            it_keys = ((key for key in x) for x in self.dic_lower if x)
            keys = chain.from_iterable(it_keys)
            keys = sorted(set(keys))
            for idx in self.indexoflist:
                for key in keys:
                    f_output.write('{}\t{}\t{}\t{}\n'.format(idx,
                                   key, repr(self.dic_top[idx][key]),
                                   repr(self.dic_lower[idx][key])))

    def _reverse_strand(self):
        curr_pos = self.record.aend-_BASES_CHECK
        fast_string = self.fasta.fetch_string(self.chrom,
                                              curr_pos, _BASES_CHECK)
        cigar_op, cigar_len = self.record.cigar[-1]
        bases = self.record.seq[-_BASES_CHECK:]

        for fast_idx in self._getindexes(fast_string):
            inverse_idx = _BASES_CHECK - fast_idx
            if (fast_idx <= _SIZE and
                    cigar_op == 0 and cigar_len >= inverse_idx and
                    bases[fast_idx:fast_idx+2] in _MINUS_STRAND_BASES):
                lst_idx = inverse_idx-2
                (self.lst_dic_lastpos[lst_idx][curr_pos+fast_idx]
                    [bases[fast_idx+1]]) += 1
                if bases[fast_idx+1] == 'A':
                    seq_for_inv = self.record.seq[-inverse_idx -
                                                  _BASE_COMP_ANALYZED:
                                                  -inverse_idx]
                    self._get_basescomp(_reverse_complement(seq_for_inv))
                    self.basecomp_lst_app(_reverse_complement(seq_for_inv))

    def _forward_strand(self):
        curr_pos = self.record.pos
        fast_string = self.fasta.fetch_string(self.chrom,
                                              curr_pos, _BASES_CHECK)
        bases = self.record.seq[:_BASES_CHECK]
        cigar_op, cigar_len = self.record.cigar[0]
        for fast_idx in self._getindexes(fast_string):
            if (fast_idx >= _SKIPBASES and
                    cigar_op == 0 and cigar_len >= fast_idx and
                    bases[fast_idx:fast_idx+2] in _PLUS_STRAND_BASES):
                (self.lst_dic_lastpos[fast_idx]
                 [curr_pos+fast_idx][bases[fast_idx]]) += 1

                if bases[fast_idx] == 'T':
                    seq = self.record.seq[fast_idx+2:
                                          _BASE_COMP_ANALYZED+fast_idx+2]
                    self._get_basescomp(seq)
                    self.basecomp_lst_app(seq)

    def write_bed(self):
        ''' every row contain chrom, start, end, top, lower, ratio'''
        fmt = '{r.chrom}\t{r.start}\t{r.end}\t\
                {r.top}\t{r.lower}\t{r.ratio}\n'.format
        with gzip.open('methyl_ms.bed.gz', 'wb') as f:
            for row in self.bed_list:
                f.write(fmt(r=row))

    def _get_basescomp(self, sequence):
        for idx, base in enumerate(sequence):
            self.basecomp_lst_dic[idx][base] += 1

    def write_basecomp(self):
        ''' dfs '''
        fmt = '{}\t{}\t{}\t{}\t{}\n'.format
        with open('methyl_ms_basecomp.txt', 'w') as f:
            f.write(fmt('pos', 'A', 'C', 'G', 'T'))
            for idx, dic in enumerate(self.basecomp_lst_dic):
                f.write(fmt(idx+1,
                            dic.get('A', 0),
                            dic.get('C', 0),
                            dic.get('G', 0),
                            dic.get('T', 0)))

    def write_basecomp_simple(self):
        with gzip.open('methyl_ms_basecomp_seqs.txt.gz', 'wb') as f:
            for line in self.basecomp_lst:
                f.write('{}\n'.format(line))


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="a bed file format with\
                        Sequences Of Interest")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...',
                        default='out_mapmethylX_bases_USER.txt')
    return parser.parse_args(argv)


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')[:3]
                chrom = input_line.pop(0).replace('chr', '')
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        yield (args.chrom, args.start, args.end)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    met_lev = Methyl_Level(args)
    last_chrom = None
    for chrom, start, end in read_bed(args):
        met_lev.reset_dict(chrom, start, end)
        last_chrom = None
        for record in samfile.fetch(chrom, start, end):
            if last_chrom and last_chrom != record.tid:
                met_lev.call_ms()
                met_lev.reset_dict(chrom, start, end)
            met_lev.update(record)
            last_chrom = record.tid
        met_lev.call_ms()
    met_lev.write_bed()
    met_lev.write_basecomp()
    met_lev.write_basecomp_simple()
    met_lev.writetofile()
    samfile.close()
    met_lev.fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
