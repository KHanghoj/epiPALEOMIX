#!/usr/bin/env python
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
from epiomix.tools.commonutils import \
    read_bed, \
    Cache

_FORW_STRAND_BASES = ['CG', 'TG']
_REV_STRAND_BASES = ['CG', 'CA']

CONV_REV = {'N': 'N',
            'C': 'G',
            'T': 'A',
            'G': 'C',
            'A': 'T'}
CONV_FORW = {'N': 'N',
             'C': 'C',
             'T': 'T',
             'G': 'G',
             'A': 'A'}


class Methyl_Level(object):
    """docstring for Methyl_Level"""
    def __init__(self, arg):
        self.arg = arg
        self._ReadBases = self.arg.ReadBases+1
        self._fasta = Cache(self.arg.FastaPath)
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.pat = re.compile('CG')
        self.rows = list()

        self.rowsapp = self.rows.append
        self.na_tup = namedtuple('row', ('pos top lower'))
        self._makeoutputfile()
        self._row_size = int(1e6)
        self._updatefunc = self._choosefunc(self.arg.Primes.lower())
        self.forw_five = self._createargdict(_FORW_STRAND_BASES, 0,
                                             self.arg.SkipFivePrime,
                                             CONV_FORW)
        self.forw_three = self._createargdict(_FORW_STRAND_BASES, 0,
                                              self.arg.SkipThreePrime,
                                              CONV_FORW)
        self.rev_five = self._createargdict(_REV_STRAND_BASES, 1,
                                            self.arg.SkipFivePrime,
                                            CONV_REV)
        self.rev_three = self._createargdict(_REV_STRAND_BASES, 1,
                                             self.arg.SkipThreePrime,
                                             CONV_REV)
        # This is all the addons
        self._count_pos_strand = [{x: 0 for x in _FORW_STRAND_BASES}
                                  for base in xrange(self._ReadBases)]
        self._count_neg_strand = [{x: 0 for x in _REV_STRAND_BASES}
                                  for base in xrange(self._ReadBases)]
        self._tempReadbase = self._ReadBases

        # in a sam/bam file everything is plus
        # strand oriented.even sequences, cigar, EVERYTHING
        # cig_idx -1 returns last cigar of the sequence from a positive strand
        # perspective. so [-1] can be the 5' end of a reverse read.

    def _createargdict(self, bases, basep, skip, convdict):
        return {'inbases': bases, 'basepos': basep,
                'skip':  skip, 'conv': convdict}

    def _choosefunc(self, prims):
        dic = {'both': self._both, 'five': self._only_5, 'three': self._only_3}
        return dic[prims]

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
        # both is only for single strand libraries as the 3' overhangs
        # are remove in double strand library construction
        if self.record.is_reverse:
            self._rightpart(**self.rev_five)
            self._leftpart(**self.rev_three)
        else:
            self._leftpart(**self.forw_five)
            self._rightpart(**self.forw_three)

    def update(self, record):
        self.record = record
        if self.record.alen < self._ReadBases:
            self._tempReadbase = self._ReadBases
            self._ReadBases = self.record.alen
        # if new read is not overlapping previous
        if self.last_end < self.record.pos:
            self._call_ms()
            self.dic_pos.clear()
            if len(self.rows) > self._row_size:
                self._writetofile()
                del self.rows[:]  # clear all data
        self._updatefunc()
        self.last_end = self.record.aend+1000  # to avoid overlapping sites.
        self._ReadBases = self._tempReadbase

    def _call_ms(self):
        for pos, basescore in sorted(self.dic_pos.iteritems()):
            # make sure overlapping genomic sites do not get counted twice
            if pos >= self.start and pos <= self.end:
                top = basescore.get('T', 0)
                low = top+basescore.get('C', 0)
                self.rowsapp(self.na_tup(pos, top, low))

    def _generate_cigar_set(self, record):
        ''' returns a set with the positions in a
        read that match/mismatch the reference '''
        readmatch = set()
        _jump_idx = 0
        cigar = record.cigar[::-1] if record.is_reverse else record.cigar
        addon = record.alen - 1 if record.is_reverse else 0
        for (cigar, count) in cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(_jump_idx, _jump_idx + count):
                    readmatch.add(abs(addon-idx))
                    if idx > self._ReadBases:
                        return readmatch
                _jump_idx += count
            elif cigar in (2, 3, 6):
                _jump_idx += count
        return readmatch

    def _prep(self, curr_pos, skip=0):
        # skip is only used in right side functions
        fast_string = self._fasta.fetch_string(self.chrom, curr_pos,
                                               self._ReadBases-skip)
        return [m.start() for m in self.pat.finditer(fast_string)]

    def _alltheread():
        ''' returns CpG hits across the entire read '''
        pass

    def _rightpart(self, inbases, basepos, skip, conv):
        curr_pos = self.record.aend-self._ReadBases
        cpg_indexes = self._prep(curr_pos, skip)
        if cpg_indexes:
            cigarmatch = self._generate_cigar_set(self.record)
            bases = self.record.seq[-self._ReadBases:]

            for fast_idx in cpg_indexes:
                inverse_idx = self._ReadBases-fast_idx
                curr_cpgcheck = bases[fast_idx:fast_idx+2]
                if ((self.record.alen-inverse_idx) in cigarmatch and
                    (self.record.alen-inverse_idx+1) in cigarmatch and
                        curr_cpgcheck in inbases):
                    # TODO:
                    # if looking at both ends of reads,
                    # we need to change the index accordingly.
                    if self.record.is_reverse:
                        (self._count_neg_strand[inverse_idx-1]
                         [curr_cpgcheck]) += 1
                    else:
                        (self._count_pos_strand[fast_idx]
                         [curr_cpgcheck]) += 1

                    (self.dic_pos[curr_pos+fast_idx+1]
                     [conv[bases[fast_idx+basepos]]]) += 1

    def _leftpart(self, inbases, basepos, skip, conv):
        curr_pos = self.record.pos
        cpg_indexes = self._prep(curr_pos)
        if cpg_indexes:
            cigarmatch = self._generate_cigar_set(self.record)
            bases = self.record.seq[:self._ReadBases]
            for fast_idx in cpg_indexes:
                inverse_idx = self._ReadBases-fast_idx
                curr_cpgcheck = bases[fast_idx:fast_idx+2]
                if (fast_idx >= skip and fast_idx in cigarmatch and
                    fast_idx+1 in cigarmatch and
                        curr_cpgcheck in inbases):

                    if self.record.is_reverse:
                        (self._count_neg_strand[inverse_idx-1]
                         [curr_cpgcheck]) += 1
                    else:
                        (self._count_pos_strand[fast_idx]
                         [curr_cpgcheck]) += 1

                    (self.dic_pos[curr_pos+fast_idx+1]
                     [conv[bases[fast_idx+basepos]]]) += 1

    def _writetofile(self):
        ''' every row contain chrom, genomicpos, top, lower, bedcoord'''
        for row in self.rows:
            self.f_output.write(self.fmt(r=row,
                                         chrom=self.chrom,
                                         bed=self.bedcoord))

    def _print_position_counts(self):
        ''' return summary for each ReadBase analyzed '''
        sys.stdout.write("Strand\tPLUS\tPLUS\tNEGA\tNEGA\n")
        sys.stdout.write("CpGPos\tTG\tCG\tCA\tCG\n")

        for base in xrange(self._ReadBases):
            curr_pos = self._count_pos_strand[base]
            curr_neg = self._count_neg_strand[base]
            sys.stdout.write('{}\t{}\t{}\t{}\t{}\n'.format(base+1,
                                                           curr_pos["TG"],
                                                           curr_pos["CG"],
                                                           curr_neg["CA"],
                                                           curr_neg["CG"]))

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
        self.fmt = "{chrom}\t{r.pos}\t{r.top}\t{r.lower}\t{bed}\n".format

    def call_final_ms(self):
        self._call_ms()
        self._writetofile()
        del self.rows[:]  # clear all data

    def reset_dict(self, chrom, start, end, bedcoord):
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = bedcoord
        del self.rows[:]  # clear all data
        self.last_end = 0

    def closefiles(self):
        self._print_position_counts()
        self.f_output.close()
        self._fasta.closefile()


def parse_helper(s):
    return s.lower()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--FastaPath', help="FastaPath: %(default)s", type=str)
    parser.add_argument('--ReadBases', help="%(default)d", type=int, default=2)
    parser.add_argument('--MinMappingQuality', help="%(default)d",
                        type=int, default=20)
    parser.add_argument('--MinAlignmentLength', help="%(default)d",
                        type=int, default=25)
    parser.add_argument('--Primes', help="%(default)s",
                        type=parse_helper, default='five',
                        choices=['both', 'five', 'three'])
    parser.add_argument('--SkipThreePrime', help="%(default)d",
                        type=int, default=0)
    parser.add_argument('--SkipFivePrime', help="%(default)d",
                        type=int, default=0)
    return parser.parse_known_args(argv)


def run(args):
    assert (args.ReadBases > args.SkipThreePrime or
            args.ReadBases > args.SkipFivePrime), \
        ("--ReadBases '{}' must be higher than both"
         " --SkipThreePrime '{}' and SkipFivePrime '{}'".format(
             args.ReadBases,
             args.SkipThreePrime,
             args.SkipFivePrime))
    samfile = pysam.AlignmentFile(args.bam, "rb")
    Met_Score = Methyl_Level(args)
    for chrom, start, end, bedcoord in read_bed(args):
        Met_Score.reset_dict(chrom, start, end, bedcoord)
        for record in samfile.fetch(chrom, start, end):
            if (record.mapq < args.MinMappingQuality or
                record.is_unmapped or
                    record.alen < args.MinAlignmentLength):
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
