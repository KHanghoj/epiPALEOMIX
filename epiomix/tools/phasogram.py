#!/usr/bin/env python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions aligning in same
orientation
'''
from __future__ import print_function
import sys
import pysam
import argparse
import gzip
from collections import defaultdict, deque
from epiomix.tools.commonutils import read_bed, \
    GC_correction


class Phasogram(GC_correction):
    """docstring for Phasogram"""
    def __init__(self, arg):
        self.arg = arg
        self.outputdic = defaultdict(int)
        self.forward_dic, self.reverse_dic = {}, {}
        GC_correction.__init__(self)

    def _call_output(self, dic):
        if self.arg.SubsetPileup > 1:
            pos_above_cutoff = [k for k, v in dic.iteritems() if v >= self.arg.SubsetPileup]
            pos_above_cutoff.sort()
        else:
            pos_above_cutoff = sorted(dic)
        sort_keys = deque(pos_above_cutoff)
        while sort_keys:
            old_pos = sort_keys.popleft()
            for curr in sort_keys:
                length = curr - old_pos
                if length > self.arg.MaxRange:
                    break
                else:
                    self.outputdic[length] += 1

    def _finddepth(self, pos, dic, corr):
        try:
            dic[pos] += corr
        except KeyError:
            dic[pos] = corr

    def update(self, record):
        if record.is_reverse:
            if record.aend <= self.end:
                corr = self._get_gc_corr_dep(record)
                self._finddepth(record.aend-1, self.reverse_dic, corr)
        else:
            if record.pos >= self.start:
                corr = self._get_gc_corr_dep(record)
                self._finddepth(record.pos, self.forward_dic, corr)

    def writetofile(self):
        fmt = '{}\t{}\n'
        with gzip.open(self.arg.outputfile, 'w') as f_output:
            for key, val in sorted(self.outputdic.iteritems()):
                f_output.write(fmt.format(key, val))
        try:
            self._fasta_dat.closefile()
        except AttributeError:
            pass

    def call(self):
        if self.forward_dic:
            self._call_output(self.forward_dic)
        if self.reverse_dic:
            self._call_output(self.reverse_dic)

    def reset(self, chrom, start, end):
        self.chrom = chrom
        self.forward_dic.clear()
        self.reverse_dic.clear()
        self.start = start-1  # 0-based as record.pos is zero based
        self.end = end


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...", type=str)
    parser.add_argument('bed', help="...", type=str)
    parser.add_argument('outputfile', help='...', type=str)
    parser.add_argument('--FastaPath', help="fastafile", type=str)
    parser.add_argument('--GCmodel', help='...', type=str, default=None)
    parser.add_argument('--SubsetPileup', help="...", type=int, default=1)
    parser.add_argument('--MaxRange', help="...", type=int, default=1500)
    parser.add_argument('--MinMappingQuality', help="...",
                        type=int, default=25)
    parser.add_argument('--MinAlignmentLength', help="...",
                        type=int, default=25)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.AlignmentFile(args.bam, "rb")
    Phaso = Phasogram(args)
    for chrom, start, end, bedcoord in read_bed(args):
        Phaso.reset(chrom, start, end)
        for record in samfile.fetch(chrom, start, end):
            if (record.mapq < args.MinMappingQuality or 
                record.is_unmapped or 
                record.is_duplicate or
                record.is_secondary or      # this is primarily for BWA MEM
                record.is_supplementary or  # this is primarily for BWA MEM
                record.is_qcfail or
                record.alen < args.MinAlignmentLength):
                continue  # do not analyze low quality records
            Phaso.update(record)
        Phaso.call()
    Phaso.writetofile()
    samfile.close()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
