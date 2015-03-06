# !/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict
from itertools import islice, izip, tee
# import gzip

_BUFFER = 0


class Cache_fasta(object):
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


class GCcorrect(object):
    """docstring for GCcorrect"""
    def __init__(self, arg):
        self.arg = arg
        self.dic_n_gc = defaultdict(lambda: defaultdict(int))
        self.dic_f_gc = defaultdict(lambda: defaultdict(int))
        self.samfile = pysam.Samfile(self.arg.bam, "rb")
        self.fasta = Cache_fasta(self.arg.fastafile)

    def getreads(self, chrom, start, end):
        self.dic_forward = {}
        self.dic_reverse = {}
        self.chrom, self.start, self.end = chrom, start, end
        records = self.samfile.fetch(chrom, start, end)
        [self._update(record.aend, self.dic_reverse) if record.is_reverse else
            self._update(record.pos, self.dic_forward) for record in records]
        if len(self.dic_forward)+len(self.dic_reverse) > 50:
            # do not use chunks with 50 or less positions.
            self._retrieve_fastaseq()
            for rl in self._read_lengths():
                for pos, rl_seq in self._short_seq(rl):
                    curr_start = pos+self.start
                    curr_end = curr_start+rl
                    gc = rl_seq.count('C')+rl_seq.count('G')
                    self.dic_n_gc[rl][gc] += 2
                    self.dic_f_gc[rl][gc] += \
                        (self.dic_forward.get(curr_start-1, 0) +
                         self.dic_reverse.get(curr_end, 0))

    def _update(self, pos, dic):
        dic[pos] = dic.get(pos, 0) + 1

    def _retrieve_fastaseq(self):
        self.length = self.end-self.start
        self.seq = self.fasta.fetch_string(self.chrom,
                                           self.start,
                                           self.length)

    def _short_seq(self, rl):
        ''' seq_fasta, length of read, return seq'''
        indices = iter(xrange(self.length - (rl - 1)))
        while True:
            try:
                idx = indices.next()
                # yield self.seq[idx:idx+rl]
                yield idx+_BUFFER, self.seq[idx+_BUFFER: idx-_BUFFER+rl]
            except StopIteration:
                break

    def _read_lengths(self):
        if self.arg.readlengths:
            with open(self.arg.readlengths, 'r') as myfile:
                for length in myfile.readlines():
                    length = int(length)
                    if length < 31:  # do not take low readlength into account
                        continue
                    else:
                        yield length
        else:
            yield self.arg.onelength

    def writetofile(self):
        ''' dfs '''
        with open(self.arg.out, 'w') as f_output:
            for rl in sorted(self.dic_n_gc):
                # for key in sorted(dic_n_gc[rl]):
                for gc in range(0, rl+1):
                    f_output.write('{}\t{}\t{}\t{}\n'.format(rl, gc,
                                   repr(self.dic_f_gc[rl][gc]),
                                   repr(self.dic_n_gc[rl][gc])))


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('mappability_bed', help="...")
    parser.add_argument('--readlengths', help="...", default=None)
    parser.add_argument('--uniqueness', help="...", type=float, default=0.9)
    parser.add_argument('--onelength', help="...", type=int, default=56)
    parser.add_argument('--out', help='...', default='out_gc_content.txt')
    return parser.parse_args(argv)


def read_bed(args):
    if args.mappability_bed:
        with open(args.mappability_bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')
                chrom = input_line.pop(0)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                score = float(input_line[-1])
                yield (chrom, start, end, score)


def main(argv):
    args = parse_args(argv)
    GC = GCcorrect(args)
    mappability = args.uniqueness
    last_chrom, last_end = '', -1
    for chrom, start, end, score in read_bed(args):
        if score >= mappability:
            if start-last_end < 0 and last_chrom == chrom:
                start = start + ((end-start)/2)
            last_chrom = chrom
            last_end = end
            GC.getreads(chrom, start, end)
    GC.writetofile()

    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
