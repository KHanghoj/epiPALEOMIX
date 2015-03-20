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
from collections import defaultdict
from itertools import chain

_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']
_BASES_CHECK = 6
_SKIPBASES = 0  # include the first base of read
_SIZE = _BASES_CHECK-_SKIPBASES-1


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
        self.lst_dic_lastpos = self._create_lst_dic(nested=True)
        self.lst_base_forward = self._create_lst_dic()
        self.dic_top = self._create_lst_dic()
        self.dic_lower = self._create_lst_dic()
        self.indexoflist = range(len(self.lst_dic_lastpos))
        self.record, self.start, self.chrom = None, None, None
        self.last_tid, self.last_pos = None, -1
        self.pat = re.compile('CG')

    def reset_dicts(self, start, chrom):
        self.lst_dic_lastpos = self._create_lst_dic(nested=True)
        self.lst_base_forward = self._create_lst_dic()
        self.last_tid, self.last_pos = None, -1
        self.chrom, self.start = chrom, start

    def _create_lst_dic(self, nested=False, size=_BASES_CHECK-1):
        if nested:
            return [defaultdict(lambda: defaultdict(int))
                    for _ in range(size)]
        else:
            return [defaultdict(int) for _ in range(size)]

    def _getindexes(self, bases_str):
        ''' returns the 0-based indeces of fasta read'''
        return [m.start() for m in self.pat.finditer(bases_str)]

    def _check_lst_dic(self):
        dic_idx = [idx for idx, x in enumerate(self.lst_dic_lastpos) if x]
        if dic_idx:
            vals = [max(self.lst_dic_lastpos[x]) for x in dic_idx]
            if vals and max(vals) < self.record.pos:
                self.ms_reverse()
                self.lst_dic_lastpos = self._create_lst_dic(nested=True)

    def update(self, record):
        self.record = record
        if self.record.tid != self.last_tid:
            # new chromosome
            if self.last_tid:
                self.ms_forward()
                self.ms_reverse()
                self.reset_dicts()

        self._check_lst_dic()

        if self.record.is_reverse:
            self._reverse_strand()
        else:
            self._forward_strand()
        self.last_tid = self.record.tid

    def call_final_ms(self):
        if self.last_pos != -1:
            self.ms_forward()
            self.ms_reverse()

    def ms_forward(self):
        self.keypos = self.last_pos-self.start
        for idx in self.indexoflist:
            tempdic_last_pos = \
                self.lst_dic_lastpos[idx].pop(self.last_pos+idx, {})
            self._call_ms(idx, self.lst_base_forward[idx],
                          tempdic_last_pos)

    def ms_reverse(self):
        for idx in self.indexoflist:
            for key in self.lst_dic_lastpos[idx].iterkeys():
                # if self.last_pos > key: ## ok because we call
                # reverse only at the end of the bed file.
                # this is the fastest solution
                # i think
                self.keypos = key-self.start
                # tempdic_last_pos = self.lst_dic_lastpos[idx].pop(key, {})
                tempdic_last_pos = self.lst_dic_lastpos[idx].get(key, {})
                self._call_ms(idx, {}, tempdic_last_pos)

    def _call_ms(self, idx, dic_forward, dic_lastpos):
        top = (dic_forward.get('T', 0) +
               dic_lastpos.get('A', 0))
        lower = (top+dic_forward.get('C', 0) +
                 dic_lastpos.get('G', 0))
        self.dic_top[idx][self.keypos+idx] += top
        self.dic_lower[idx][self.keypos+idx] += lower

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
        fast_idx = self._getindexes(fast_string)
        if fast_idx:
            bases = self.record.seq[-_BASES_CHECK:]
            read_idx = [x for x in fast_idx if
                        bases[x:x+2] in _MINUS_STRAND_BASES and
                        x < _SIZE]  # x <= _SIZE]
            if read_idx:
                max_pos = _BASES_CHECK - min(read_idx)
                cigar_op, cigar_len = self.record.cigar[-1]
                if (cigar_op == 0) and (cigar_len >= max_pos):
                    for base_idx in read_idx:
                        (self.lst_dic_lastpos[_BASES_CHECK-base_idx-2]
                            [curr_pos+base_idx+1]
                            [bases[base_idx+1]]) += 1

    def _forward_strand(self):
        curr_pos = self.record.pos
        fast_string = self.fasta.fetch_string(self.chrom,
                                              curr_pos, _BASES_CHECK)
        fast_idx = self._getindexes(fast_string)
        if fast_idx:
            bases = self.record.seq[:_BASES_CHECK]
            read_idx = [x for x in fast_idx if
                        bases[x:x+2] in _PLUS_STRAND_BASES and
                        x >= _SKIPBASES]
            if read_idx:
                max_pos = max(read_idx)+2
                # because we look at the two following bases.
                cigar_op, cigar_len = self.record.cigar[0]
                if cigar_op == 0 and cigar_len >= max_pos:
                    if curr_pos != self.last_pos and self.last_pos != -1:
                        self.ms_forward()
                        # RESET DICT:
                        self.lst_base_forward = self._create_lst_dic()

                    for base_idx in read_idx:
                        self.lst_base_forward[base_idx][bases[base_idx]] += 1
                    self.last_pos = curr_pos


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
    for chrom, start, end in read_bed(args):
        met_lev.reset_dicts(start, chrom)
        for record in samfile.fetch(chrom, start, end):
            met_lev.update(record)
        met_lev.call_final_ms()
    met_lev.writetofile()
    samfile.close()
    met_lev.fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
