# !/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict
from epiomix_commonutils import \
    strtobool, Cache, corr_fasta_chr
_BUFFER = 2


class GCcorrect(object):
    """docstring for GCcorrect"""
    def __init__(self, arg):
        self.arg = arg
        self.reads_gc, self.reference_gc = defaultdict(int), defaultdict(int)
        self.samfile = pysam.Samfile(self.arg.BamPath, "rb")
        self.fasta = Cache(self.arg.FastaPath)
        self.rl = self.arg.ReadLength

    def getreads(self, chrom, start, end):
        self.reads_forw, self.reads_back = {}, {}
        records = self.samfile.fetch(chrom, start, end)
        # they need to be 0-based for fetching fasta seq:
        self.chrom = corr_fasta_chr(self.arg, chrom)
        self.start, self.end = start-1, end-1
        for record in records:
            if record.is_reverse:
                self._update(record.aend-1, self.reads_back)
            else:
                self._update(record.pos, self.reads_forw)

        if self.reads_forw and self.reads_back:
            for rela_pos, gc in self._short_seq(self.rl):
                curr_start = rela_pos+self.start
                curr_end = curr_start+self.rl
                self.reference_gc[gc] += 2
                self.reads_gc[gc] += (self.reads_forw.pop(curr_start, 0) +
                                      self.reads_back.pop(curr_end, 0))

    def _update(self, pos, dic):
        try:
            dic[pos] += 1
        except KeyError:
            dic[pos] = 1

    def _retrieve_fastaseq(self):
        length = self.end-self.start
        seq = self.fasta.fetch_string(self.chrom, self.start, length)
        return length, seq

    def _short_seq(self, rl):
        ''' seq_fasta, length of read, return seq'''
        region_size, region_seq = self._retrieve_fastaseq()
        for idx in xrange(region_size - (rl - 1)):
            seq = region_seq[idx+_BUFFER: idx+rl-_BUFFER]
            yield (idx+_BUFFER, (seq.count('C')+seq.count('G')))

    def writetofile(self):
        ''' dfs '''
        # out = sys.stdout.write
        fmt = '{}\t{}\t{}\t{}\n'.format
        with open(self.arg.OutputFile, 'w') as f:
            for gc in range(0, self.rl+1):
                # out(fmt(str(self.rl), str(gc), str(self.reads_gc[gc]),
                #     str(self.reference_gc[gc])))
                f.write(fmt(str(self.rl), str(gc), str(self.reads_gc[gc]),
                        str(self.reference_gc[gc])))
        self.fasta.closefile()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='GCcorrection')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('OutputFile', type=str)
    parser.add_argument('--FastaPath', type=str)
    parser.add_argument('--MappabilityPath', type=str)
    parser.add_argument('--ReadLength', help="...", type=int)
    parser.add_argument('--MappaUniqueness', help="...", type=float)
    parser.add_argument('--FastaChromType')
    parser.add_argument('--BamChromType')
    parser.add_argument('--OffSet', type=str,
                        help='the offsetfile found by the midnode')
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    return parser.parse_known_args(argv)


def read_mappa_W(args):
    with open(args.MappabilityPath, 'r') as mappafile:
        for line in mappafile:
            input_line = (line.rstrip('\n')).split('\t')
            chrom, start, end = input_line[:3]
            if 'chr' not in chrom:
                chrom = 'chr{}'.format(chrom)
            score = float(input_line[-1])
            yield (str(chrom), int(start), int(end), score)


def read_mappa_WO(args):
    with open(args.MappabilityPath, 'r') as mappafile:
        for line in mappafile:
            input_line = (line.rstrip('\n')).split('\t')
            chrom, start, end = input_line[:3]
            if 'chr' in chrom:
                chrom = chrom.replace('chr', '')
            score = float(input_line[-1])
            yield (str(chrom), int(start), int(end), score)


def run(args):
    read_bed = read_mappa_W if strtobool(args.BamChromType) else read_mappa_WO
    args.FastaChromType = strtobool(args.FastaChromType)
    if args.OffSet:
        args.ReadLength += int(open(args.OffSet, 'r').read().strip())
    GC = GCcorrect(args)
    mappability = args.MappaUniqueness
    last_chrom, last_end = '', -1
    for chrom, start, end, score in read_bed(args):
        if score >= mappability:  # and chrom in 'chr1' or chrom in '1':
            # because chunks can overlap with 50%
            if start-last_end < 0 and last_chrom == chrom:
                start = start + ((end-start)/2)
            last_chrom, last_end = chrom, end
            GC.getreads(chrom, start, end)
    GC.writetofile()


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
