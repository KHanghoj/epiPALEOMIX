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
        self.dic_n_gc = defaultdict(int)
        self.dic_f_gc = defaultdict(int)
        self.samfile = pysam.Samfile(self.arg.BamPath, "rb")
        self.fasta = Cache(self.arg.FastaPath)
        self.rl = self.arg.ReadLength

    def getreads(self, chrom, start, end):
        self.dic_forward, self.dic_reverse = {}, {}
        records = self.samfile.fetch(chrom, start, end)
        # they need to be 0-based for fetching fasta seq:
        self.chrom = corr_fasta_chr(self.arg, chrom)
        self.start, self.end = start-1, end-1
        [self._update(record.aend-1, self.dic_reverse)
            if record.is_reverse else
            self._update(record.pos, self.dic_forward) for record in records]
        if self.dic_forward and self.dic_reverse:
            # do not use chunks with 50 or less positions.
            for rela_pos, gc in self._short_seq(self.rl):
                curr_start = rela_pos+self.start
                curr_end = curr_start+self.rl
                self.dic_n_gc[gc] += 2
                self.dic_f_gc[gc] += (self.dic_forward.pop(curr_start, 0) +
                                      self.dic_reverse.pop(curr_end, 0))

    def _update(self, pos, dic):
        try:
            dic[pos] += 1
        except KeyError:
            dic[pos] = 1

    def _retrieve_fastaseq(self):
        length = self.end-self.start
        seq = self.fasta.fetch_string(self.chrom, self.start-1, length)
        return(length, seq)

    def _short_seq(self, rl):
        ''' seq_fasta, length of read, return seq'''
        region_size, region_seq = self._retrieve_fastaseq()
        indices = iter(xrange(region_size - (rl - 1)))
        while True:
            try:
                idx = indices.next()
                tempseq = region_seq[idx+_BUFFER: idx-_BUFFER+rl]
                yield (idx+_BUFFER, (tempseq.count('C')+tempseq.count('G')))
            except StopIteration:
                break

    def writetofile(self):
        ''' dfs '''
        out = sys.stdout.write
        fmt = '{}\t{}\t{}\t{}\n'.format
        for gc in range(0, self.rl+1):
            out(fmt(str(self.rl), str(gc), str(self.dic_f_gc[gc]),
                str(self.dic_n_gc[gc])))


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='GCcorrection')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('--FastaPath', type=str)
    parser.add_argument('--MappabilityPath', type=str)
    parser.add_argument('--ReadLength', help="...", type=int)
    parser.add_argument('--MappaUniqueness', help="...", type=float)
    parser.add_argument('--FastaChromType', dest='FastaChromType')
    parser.add_argument('--BamChromType', dest='BamChromType')
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    return parser.parse_known_args(argv)


def read_mappa_W(args):
    with open(args.MappabilityPath, 'r') as myfile:
        for line in myfile:
            input_line = (line.rstrip('\n')).split('\t')
            chrom, start, end = input_line[:3]
            if 'chr' not in chrom:
                chrom = 'chr{}'.format(chrom)
            score = float(input_line[-1])
            yield (str(chrom), int(start), int(end), score)


def read_mappa_WO(args):
    with open(args.MappabilityPath, 'r') as myfile:
        for line in myfile:
            input_line = (line.rstrip('\n')).split('\t')
            chrom = input_line.pop(0).replace('chr', '')
            if 'chr' in chrom:
                chrom = chrom.replace('chr', '')
            start = int(input_line.pop(0))
            end = int(input_line.pop(0))
            score = float(input_line[-1])
            yield (chrom, start, end, score)


def run(args):
    read_bed = read_mappa_W if strtobool(args.BamChromType) else read_mappa_WO
    args.FastaChromType = strtobool(args.FastaChromType)
    GC = GCcorrect(args)
    mappability = args.MappaUniqueness
    last_chrom, last_end = '', -1
    for chrom, start, end, score in read_bed(args):
        if score >= mappability:
            # because chunks can overlap
            if start-last_end < 0 and last_chrom == chrom:
                start = start + ((end-start)/2)
            last_chrom = chrom
            last_end = end
            GC.getreads(chrom, start, end)
    GC.writetofile()


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))