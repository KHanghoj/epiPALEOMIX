#!/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
import gzip
from collections import deque
# from itertools import islice, izip, tee
from os import remove
_MAXLEN = 500


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


class Write_Depth(object):
    """docstring for Write_Depth"""
    def __init__(self, arg, seq_len=_MAXLEN):
        self.arg = arg
        self._fasta = Cache(self.arg.fastafile)
        self._f_output, self._model, self._GC_model_len = None, None, None
        self._seq_len = int(seq_len)
        self._last_pos = -self._seq_len
        self._zeros = [0] * self._seq_len
        self._corrected_depth = deque(self._zeros, maxlen=self._seq_len)
        self._genomic_positions = deque(maxlen=self._seq_len)
        self._cache_list = list()
        self._output_size = self._seq_len * 10
        self.counter, self.deque_idx = 0, 0
        self._output_fmt = '{0}\t{1}\t{2}\n'
        self._makeoutputfile()
        self._gcmodel_ini()

    def _gcmodel_ini(self):
        if self.arg.gcmodel:
            with open(self.arg.gcmodel, 'r') as f:
                self._model = [float(line.rstrip('\n').split('\t')[-1])
                               for line in f]
                self._GC_model_len = len(self._model)
        else:
            self._GC_model_len = 0  # this is default if not assign
            self._model = [1]*(self._GC_model_len+1)

    def _get_gc_corr_dep(self, pos):
        fasta_str = self._fasta.fetch_string(self._present_chrom,
                                             pos, self._GC_model_len)
        gc_idx = fasta_str.count('G')+fasta_str.count('C')
        return self._model[gc_idx]

    def update_depth(self, record, start, chrom):
        self._present_chrom, self.start = chrom, start
        curr_pos = record.pos
        self.jump = (record.pos - self._last_pos)
        if self.jump > 0:
            if self.jump >= self._seq_len:
                self.jump = self._seq_len
                if self._genomic_positions:
                    self._append_cache_lst()
                self._corrected_depth = deque(self._zeros,
                                              maxlen=self._seq_len)
                self._genomic_positions = deque(((pos+1) for pos in
                                                xrange(curr_pos,
                                                       curr_pos +
                                                       self._seq_len)),
                                                maxlen=self._seq_len)
            else:
                self._append_cache_lst()
                self._corrected_depth.extend([0]*self.jump)
                position_values = curr_pos+self._seq_len
                self._genomic_positions.extend((pos+1) for pos in
                                               xrange(position_values,
                                                      position_values +
                                                      self.jump))

        if record.is_reverse:
            corr_depth = self._get_gc_corr_dep(record.aend-self._GC_model_len)
        else:
            corr_depth = self._get_gc_corr_dep(record.pos)

        self.deque_idx = 0
        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(self.deque_idx, self.deque_idx + count):
                    self._corrected_depth[idx] += corr_depth
                self.deque_idx += count
            elif cigar in (2, 3, 6):
                self.deque_idx += count
        self._last_pos = record.pos

    def adjust_deques(self, jump):
        while jump:
            yield (self._genomic_positions.popleft(),
                   self._corrected_depth.popleft())
            jump -= 1

    def _append_cache_lst(self):
        for pos, dep in self.adjust_deques(self.jump):
            if dep and pos >= self.start:
                # self.f_output.write(self._output_fmt.format(pos, dep))
                self._cache_list.append((pos, dep))
                self.counter += 1
        if self.counter > self._output_size:
            self._write_to_file()

    def _write_to_file(self):
        ''' dfs '''
        if self._cache_list:
            for p, d in iter(self._cache_list):
                self.f_output.write(self._output_fmt.
                                    format(self._present_chrom, p, d))
            self._cache_list = list()
            self.counter = 0

    def call_final(self):
        self._append_cache_lst()
        self._write_to_file()

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        try:
            remove(self.arg.out)
            self.f_output = gzip.open(self.arg.out, 'ab')
        except OSError:
            self.f_output = gzip.open(self.arg.out, 'ab')

    def reset_deques(self):
        self._corrected_depth.clear()
        self._genomic_positions.clear()
        self._last_pos = -self._seq_len

    def closefile(self):
        self.f_output.close()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('--gcmodel', help="fastafile")
    parser.add_argument('--bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", default=None)
    parser.add_argument('--end', help="...", default=None)
    parser.add_argument('--out', help='...', default='out_corr_depth.txt.gz')
    return parser.parse_args(argv)


def read_bed(args, chromtype=''):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = line.rstrip('\n').split('\t')
                chrom = input_line.pop(0).replace('chr', chromtype)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        yield (args.chrom, int(args.start), int(args.end))


def main(argv):
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    corrected_depth = Write_Depth(args)
    for chrom, start, end in read_bed(args):
        last_tid = None
        corrected_depth.reset_deques()
        for record in samfile.fetch(chrom, start, end):
            if record.tid != last_tid:
                if last_tid:
                    corrected_depth.call_final()
                    corrected_depth.reset_deques()
            corrected_depth.update_depth(record, start, chrom)
            last_tid = record.tid
        corrected_depth.call_final()
    corrected_depth.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
