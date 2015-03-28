# !/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
import tempfile
from collections import defaultdict
from itertools import repeat, izip
import os
import shutil
import multiprocessing as mp
import subprocess as sp

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
        self.dic_n_gc = defaultdict(int)
        self.dic_f_gc = defaultdict(int)
        self.samfile = pysam.Samfile(self.arg.bam, "rb")
        self.fasta = Cache_fasta(self.arg.fastafile)

    def getreads_old(self, chrom, start, end):
        self.dic_forward = {}
        self.dic_reverse = {}
        records = self.samfile.fetch(chrom, start, end)
        # they need to be 0-based for fetching fasta seq:
        self.chrom, self.start, self.end = chrom, start-1, end-1

        [self._update(record.aend-1, self.dic_reverse)
            if record.is_reverse else
            self._update(record.pos, self.dic_forward) for record in records]
        if len(self.dic_forward)+len(self.dic_reverse) > 50:
            # do not use chunks with 50 or less positions.
            for rl in self._read_lengths():
                for rela_pos, gc in self._short_seq(rl):
                    curr_start = rela_pos+self.start
                    curr_end = curr_start+rl
                    self.dic_n_gc[rl][gc] += 2
                    self.dic_f_gc[rl][gc] += \
                        (self.dic_forward.pop(curr_start, 0) +
                         self.dic_reverse.pop(curr_end, 0))

    def getreads(self, chrom, start, end, rl):
        self.dic_forward, self.dic_reverse = {}, {}
        self.rl = rl
        records = self.samfile.fetch(chrom, start, end)
        # they need to be 0-based for fetching fasta seq:
        self.chrom, self.start, self.end = chrom, start-1, end-1
        [self._update(record.aend-1, self.dic_reverse)
            if record.is_reverse else
            self._update(record.pos, self.dic_forward) for record in records]
        if len(self.dic_forward)+len(self.dic_reverse) > 50:
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
        seq = self.fasta.fetch_string(self.chrom,
                                      self.start-1,
                                      length)
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

    # def _read_lengths(self):
    #    # do not take low readlength into account
    #     if self.arg.readlengths:
    #         with open(self.arg.readlengths, 'r') as myfile:
    #             for length in myfile.readlines():
    #                 length = int(length)
    #                 if length < 31:
    #                     continue
    #                 else:
    #                     yield length
    #     else:
    #         yield self.arg.onelength

    def writetofile(self, tempf):
        ''' dfs '''
        pref = os.path.basename(self.arg.bam)
        f_output = tempfile.NamedTemporaryFile(prefix=pref,
                                               suffix='_GCCORRECT',
                                               delete=False,
                                               dir=tempf)
        for gc in range(0, self.rl+1):
            f_output.write('{}\t{}\t{}\t{}\n'.format(self.rl, gc,
                           repr(self.dic_f_gc[gc]),
                           repr(self.dic_n_gc[gc])))
        f_output.close()


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
            for line in myfile:
                input_line = (line.rstrip('\n')).split('\t')
                chrom = input_line.pop(0)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                score = float(input_line[-1])
                yield (chrom, start, end, score)


def _run(args, rl, tempf):
    GC = GCcorrect(args)
    mappability = args.uniqueness
    last_chrom, last_end = '', -1
    for chrom, start, end, score in read_bed(args):
        if score >= mappability:
            # because chunks can overlap
            if start-last_end < 0 and last_chrom == chrom:
                start = start + ((end-start)/2)
            last_chrom = chrom
            last_end = end
            GC.getreads(chrom, start, end, rl)
    GC.writetofile(tempf)


def run_helper(argum):
    # args, rl, tempf = argum
    # return _run(args, rl, tempf)
    return _run(*argum)


def main(argv):
    tempf = tempfile.mkdtemp(dir=os.getcwd())
    pool = mp.Pool(processes=2)
    args = parse_args(argv)
    print(args)
    rl = [31, 56]
    # repeat holds the two string argument together
    # enables me two use the *argum to make it more
    # versatile
    params = izip(repeat(args), rl, repeat(tempf))
    pool.map_async(run_helper, params)
    pool.close()
    pool.join()
    # Rscript, scriptname, inputfolder, pattern
    cmd_str = ['Rscript', 'plot_gc.R', tempf, os.path.basename(args.bam)]
    sp.check_call(cmd_str)
    shutil.rmtree(tempf)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
