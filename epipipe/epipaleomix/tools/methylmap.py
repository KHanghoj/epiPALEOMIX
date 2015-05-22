# !/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
'''
from __future__ import print_function
import sys
import pysam
import argparse
import re
import gzip
from os.path import exists, splitext
from shutil import move
from collections import defaultdict, namedtuple
from epipaleomix.tools.commonutils import \
    read_bed, \
    Cache, \
    corr_chrom
_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']
SSorDS = {
    'SS': {'inbases': _PLUS_STRAND_BASES, 'basepos': 0},
    'DS': {'inbases': _MINUS_STRAND_BASES, 'basepos': 1}
}


class Methyl_Level(object):
    """docstring for Methyl_Level"""
    def __init__(self, arg):
        self.arg = arg
        self._ReadBases = self.arg.ReadBases
        self._skip_five = self.arg.SkipFivePrime
        self._skip_three = self.arg.SkipThreePrime
        self._fasta = Cache(self.arg.FastaPath)
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.pat = re.compile('CG')
        self.rows = list()
        self.rowsapp = self.rows.append
        self.na_tup = namedtuple('row', ('pos top lower'))
        self._makeoutputfile()
        self._row_size = int(1e6)
        self._init_dics()
        self._updatefunc = self._choosefunc()

    def _init_dics(self):
        lib_info = SSorDS[self.arg.LibraryConstruction]
        forwdic = {'inbases': _PLUS_STRAND_BASES, 'basepos': 0}
        self.forw_five = merge_dics(forwdic, {'skip': self._skip_five})
        self.forw_three =  merge_dics(forwdic, {'skip':self._skip_three})
        self.rev_five = merge_dics(lib_info, {'skip': self._skip_five})
        self.rev_three = merge_dics(lib_info, {'skip': self._skip_three})
        # in a sam/bam file everything is plus strand oriented.even sequences, cigar, EVERYTHING
        # cig_idx -1 returns last cigar of the sequence from a positive strand
        # perspective. so [-1] can be the 5' end of a reverse read.

    def _choosefunc(self):
        dic={'both': self._both, 'five': self._only_5, 'three': self._only_3}
        return dic[self.arg.Primes.lower()]
     
    def _only_5(self):
        if self.record.is_reverse:
            self._rightpart(**self.rev_five)
        else:
            self._leftpart(**self.forw_five)

    def _only_3(self):
        if self.record.is_reverse:
            self._leftpart(**self.rev_three)
        else:
            self._rightpart(**self.forw_three)

    def _both(self):
        # both is only for single strand libs
        # this is because the 3' overhangs are remove in double strand library
        # preparation.
        if self.record.is_reverse:
            self._rightpart(**self.rev_five)
            self._leftpart(**self.rev_three)
        else:
            self._leftpart(**self.forw_five)
            self._rightpart(**self.forw_three)

    def reset_dict(self, chrom, start, end, bedcoord):
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = bedcoord
        del self.rows[:]  # clear all data
        self.last_end = 0

    def _getindexes(self, bases_str):
        ''' returns the 0-based indeces of fasta read'''
        return (m.start() for m in self.pat.finditer(bases_str))

    def update(self, record):
        self.record = record
        if self.last_end < self.record.pos:
            self._call_ms()
            self.dic_pos.clear()
            if len(self.rows) > self._row_size:
                self._writetofile()
        self._updatefunc()
        self.last_end = self.record.aend

    def _call_ms(self):
        for pos, basescore in sorted(self.dic_pos.iteritems()):
            if pos >= self.start and pos <= self.end:  # make sure overlapping genomic sites do not get counted twice
                top = basescore.get('T', 0)+basescore.get('A', 0)
                low = top+basescore.get('C', 0)+basescore.get('G', 0)
                self.rowsapp(self.na_tup(pos, top, low))

    def call_final_ms(self):
        self._call_ms()
        self._writetofile()

    def _prep(self, curr_pos, skip=0):
        ## skip is only used in right side functions
        # skip is 1-based always
        fast_string = self._fasta.fetch_string(self.chrom, curr_pos,
                                               self._ReadBases-skip)
        for fast_idx in self._getindexes(fast_string):
            yield fast_idx, self._ReadBases - fast_idx

    def _rightpart(self, inbases, basepos, skip):
        curr_pos = self.record.aend-self._ReadBases
        cigar_op, cigar_len = self.record.cigar[-1]
        bases = self.record.seq[-self._ReadBases:]
        for fast_idx, inverse_idx in self._prep(curr_pos, skip):
            if (cigar_op == 0 and cigar_len >= inverse_idx and
                    bases[fast_idx:fast_idx+2] in inbases):
                self.dic_pos[curr_pos+fast_idx+1][bases[fast_idx+basepos]] += 1
        
    def _leftpart(self, inbases, basepos, skip):
        curr_pos = self.record.pos
        cigar_op, cigar_len = self.record.cigar[0]
        bases = self.record.seq[:self._ReadBases]
        for fast_idx, _ in self._prep(curr_pos):
            if (fast_idx >= skip and cigar_op == 0 and
                    cigar_len >= fast_idx and
                    bases[fast_idx:fast_idx+2] in inbases):
                self.dic_pos[curr_pos+fast_idx+1][bases[fast_idx+basepos]] += 1            

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
        self.fmt = "{chr}\t{r.pos}\t{r.top}\t{r.lower}\t{bed}\n".format

    def _writetofile(self):
        ''' every row contain chrom, genomicpos, top, lower, bedcoord'''
        for row in self.rows:
            self.f_output.write(self.fmt(r=row, chr=self.chrom,
                                         bed=self.bedcoord))
        del self.rows[:]  # EMPTY LIST AFTER WRITING TO FILE

    def closefiles(self):
        try:
            self._fasta.closefile()
        except AttributeError:
            pass
        self.f_output.close()


def merge_dics(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--FastaPrefix', dest='FastaPrefix')
    parser.add_argument('--BamPrefix', dest='BamPrefix')
    parser.add_argument('--ReadBases', help="..", type=int, default=6)
    parser.add_argument('--SkipBases', help="..", type=int, default=0)
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    parser.add_argument('--LibraryConstruction', help="..", type=str,
                        choices=['DS', 'SS'])
    parser.add_argument('--Primes', help="..", type=str,
                        choices=['both', 'five', 'three'])
    parser.add_argument('--SkipThreePrime', help="..", type=int, default=0)
    parser.add_argument('--SkipFivePrime', help="..", type=int, default=0)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.Samfile(args.bam, "rb")
    Met_Score = Methyl_Level(args)
    for chrom, start, end, bedcoord in read_bed(args):
        Met_Score.reset_dict(corr_chrom(args.FastaPrefix, chrom),
                             start, end, bedcoord)
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < args.MinMappingQuality:
                continue  # do not analyze low quality records
            Met_Score.update(record)
        Met_Score.call_final_ms()
    Met_Score.closefiles()
    samfile.close()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
