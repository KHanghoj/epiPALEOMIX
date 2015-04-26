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
from epiomix_commonutils import read_bed_W, \
    read_bed_WO, strtobool, Cache, corr_fasta_chr
_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']


class Methyl_Level(object):
    """docstring for Methyl_Level"""
    def __init__(self, arg):
        self.arg = arg
        self._ReadBases = self.arg.ReadBases
        self._SkipBases = self.arg.SkipBases
        self._Size = self._ReadBases-self._SkipBases-1
        self._fasta = Cache(self.arg.FastaPath)
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.pat = re.compile('CG')
        self.rows = list()
        self.rowsapp = self.rows.append
        self.na_tup = namedtuple('row', ('pos top lower'))
        self._makeoutputfile()

    def reset_dict(self, chrom, start, end):
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = '{}_{}_{}'.format(self.chrom, self.start, self.end)
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
            if len(self.rows) > int(1e6):
                self._writetofile()

        if self.record.is_reverse:
            self._minus_or_threeprime()
        else:
            self._forward_strand()
        self.last_end = record.aend

    def update_ss(self, record):
        self.record = record
        if self.last_end < self.record.pos:
            self._call_ms()
            self.dic_pos.clear()
            if len(self.rows) > int(1e6):
                self._writetofile()
        self._minus_or_threeprime(inbases=_PLUS_STRAND_BASES, libtype=0)
        self._forward_strand()
        self.last_end = record.aend

    def _call_ms(self):
        for pos, basescore in sorted(self.dic_pos.iteritems()):
            top = basescore.get('T', 0)+basescore.get('A', 0)
            low = top+basescore.get('C', 0)+basescore.get('G', 0)
            self.rowsapp(self.na_tup(pos+1, top, low))

    def call_final_ms(self):
        self._call_ms()
        self._writetofile()

    # def _reverse_strand(self):
    #     curr_pos = self.record.aend-self._ReadBases
    #     fast_string = self._fasta.fetch_string(self.chrom,
    #                                            curr_pos, self._ReadBases)
    #     cigar_op, cigar_len = self.record.cigar[-1]
    #     bases = self.record.seq[-self._ReadBases:]
    #     for fast_idx in self._getindexes(fast_string):
    #         inverse_idx = self._ReadBases - fast_idx
    #         if (fast_idx <= self._Size and
    #                 cigar_op == 0 and cigar_len >= inverse_idx and
    #                 bases[fast_idx:fast_idx+2] in _MINUS_STRAND_BASES):
    #             self.dic_pos[curr_pos+fast_idx][bases[fast_idx+1]] += 1

    def _minus_or_threeprime(self, inbases=_MINUS_STRAND_BASES, libtype=1):
        curr_pos = self.record.aend-self._ReadBases
        fast_string = \
            self._fasta.fetch_string(self.chrom, curr_pos, self._ReadBases)
        cigar_op, cigar_len = self.record.cigar[-1]
        bases = self.record.seq[-self._ReadBases:]
        for fast_idx in self._getindexes(fast_string):
            inverse_idx = self._ReadBases - fast_idx
            if (fast_idx <= self._Size and cigar_op == 0 and
                    cigar_len >= inverse_idx and
                    bases[fast_idx:fast_idx+2] in inbases):
                self.dic_pos[curr_pos+fast_idx][bases[fast_idx+libtype]] += 1

    def _forward_strand(self):
        curr_pos = self.record.pos
        fast_string = \
            self._fasta.fetch_string(self.chrom, curr_pos, self._ReadBases)
        bases = self.record.seq[:self._ReadBases]
        cigar_op, cigar_len = self.record.cigar[0]
        for fast_idx in self._getindexes(fast_string):
            if (fast_idx >= self._SkipBases and
                    cigar_op == 0 and cigar_len >= fast_idx and
                    bases[fast_idx:fast_idx+2] in _PLUS_STRAND_BASES):
                self.dic_pos[curr_pos+fast_idx][bases[fast_idx]] += 1

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
        # header = '#chrom\tgenomicpos\tdeaminated\ttotal\tbedcoord\n'
        # self.f_output.write(header)
        self.fmt = "{chr}\t{r.pos}\t{r.top}\t{r.lower}\t{bed}\n".format

    def _writetofile(self):
        ''' every row contain chrom, genomicpos, top, lower, bedcoord'''
        for row in self.rows:
            self.f_output.write(self.fmt(r=row, chr=self.chrom,
                                         bed=self.bedcoord))
        del self.rows[:]  # EMPTY LIST AFTER WRITING TO FILE

    def closefiles(self):
        self._fasta.closefile()
        self.f_output.close()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--FastaChromType', dest='FastaChromType')
    parser.add_argument('--BamChromType', dest='BamChromType')
    parser.add_argument('--ReadBases', help="..", type=int, default=6)
    parser.add_argument('--SkipBases', help="..", type=int, default=0)
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    parser.add_argument('--LibraryConstruction', help="..", type=str,
                        choices=['DS', 'SS'])
    return parser.parse_known_args(argv)


def run(args):
    read_bed = read_bed_W if strtobool(args.BamChromType) else read_bed_WO
    args.FastaChromType = strtobool(args.FastaChromType)
    samfile = pysam.Samfile(args.bam, "rb")
    Met_Score = Methyl_Level(args)
    update = (Met_Score.update if args.LibraryConstruction == 'DS'
              else Met_Score.update_ss)
    for chrom, start, end in read_bed(args):
        Met_Score.reset_dict(corr_fasta_chr(args, chrom), start, end)
        for record in samfile.fetch(chrom, start, end):
            update(record)
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
