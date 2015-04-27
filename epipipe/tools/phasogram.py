#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
import gzip
from collections import defaultdict, deque
from epipipe.tools.commonutils import read_bed_W, \
    read_bed_WO, strtobool, GC_correction, corr_fasta_chr # , Cache


class Phasogram(GC_correction):
    """docstring for Phasogram"""
    def __init__(self, arg):
        self.arg = arg
#         self._fasta = Cache(self.arg.FastaPath)
        GC_correction.__init__(self)
        self.outputdic = defaultdict(int)
        self.forward_dic, self.reverse_dic = {}, {}
        self._GC_model_len = 0
        self._GCmodel_ini()

    def _call_output(self, dic, max_range=None, max_size=None):
        if dic:
            max_range = self.arg.MaxRange if max_range is None else max_range
            max_size = self.arg.MaxRange if max_size is None else max_range
            sort_keys = deque(sorted(dic.iterkeys()))
            max_key = sort_keys[-1]
            while max_key - sort_keys[0] > max_range:
                old_pos = sort_keys.popleft()
                old_count = dic.pop(old_pos, 0)
                for current in sort_keys:
                    length = current - old_pos
                    if length >= max_size:
                        break  # break the for loop as current pos is too long.
                    if old_count >= self.arg.SubsetPileup:
                        self.outputdic[length] += 1
                        
    def _finddepth(self, pos, dic, gc_pos):
        corr = self._get_gc_corr_dep(gc_pos)
        try:
            dic[pos] += corr
        except KeyError:
            dic[pos] = corr
            self._call_output(dic)

    
    def update(self, record):
        if record.is_reverse:
            self._finddepth(record.aend, self.reverse_dic, record.aend-self._GC_model_len)
        else:
            self._finddepth(record.pos, self.forward_dic, record.pos)

    def writetofile(self):
        with gzip.open(self.arg.outputfile, 'w') as f_output:
            for key in sorted(self.outputdic.iterkeys()):
                f_output.write('{}\t{}\n'.format(key, self.outputdic[key]))
        try:
            self._fasta_dat.closefile()
        except AttributeError:
            pass

    def call(self):
        self._call_output(self.forward_dic, max_range=0, max_size=0)
        self._call_output(self.reverse_dic, max_range=0, max_size=0)

    def reset(self, chrom):
        self.chrom = chrom
        self.forward_dic, self.reverse_dic = {}, {}


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...", type=str)
    parser.add_argument('bed', help="...", type=str)
    parser.add_argument('outputfile', help='...', type=str)
    parser.add_argument('--FastaPath', help="fastafile", type=str)
    parser.add_argument('--GCmodel', help='...', type=str, default=None)
    parser.add_argument('--SubsetPileup', help="...", type=int, default=3)
    parser.add_argument('--MaxRange', help="...", type=int, default=1000)
    parser.add_argument('--FastaChromType', help="...")
    parser.add_argument('--BamChromType', help="...")
    parser.add_argument('--MinMappingQuality', help="...", type=int,
                        default=25)
    return parser.parse_known_args(argv)


def run(args):
    read_bed = read_bed_W if strtobool(args.BamChromType) else read_bed_WO
    args.FastaChromType = strtobool(args.FastaChromType) if args.GCmodel else True
#    args.FastaChromType = strtobool(args.FastaChromType)
    samfile = pysam.Samfile(args.bam, "rb")
    Phaso = Phasogram(args)
    for chrom, start, end in read_bed(args):
        Phaso.reset(corr_fasta_chr(args, chrom))
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < args.MinMappingQuality:
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
