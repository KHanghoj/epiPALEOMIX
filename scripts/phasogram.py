#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict, deque

_MAX_SIZE = 3000
_MINMAPQUALI = 25
_MIN_COVERAGE = 3


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


class Phasogram(object):
    """docstring for Phasogram"""
    def __init__(self, arg):
        self.arg = arg
        self._fasta = Cache(self.arg.fastafile)
        self.outputdic = defaultdict(int)
        self.forward_dic = {}
        self.reverse_dic = {}
        self.last_tid = None
        self._GC_model_len = None
        self._gcmodel_ini()

    def _call_output(self, dic, max_lst_range=_MAX_SIZE, max_size=_MAX_SIZE):
        if dic:
            sort_keys = deque(sorted(dic.iterkeys()))
            mdic = sort_keys[-1]
            while mdic - sort_keys[0] > max_lst_range:
                old_pos = sort_keys.popleft()
                old_count = dic.pop(old_pos, None)
                for current in sort_keys:
                    length = current - old_pos
                    if length >= max_size:
                        break
                        # do no know if break or contiune
                    if old_count >= _MIN_COVERAGE:
                        self.outputdic[length] += 1

    def update(self, record, chrom):
        self.chrom = chrom
        if record.is_reverse:
            curr_pos, temp_dic = record.aend, self.reverse_dic
            corr_depth = self._get_gc_corr_dep(record.aend-self._GC_model_len)
        else:
            curr_pos, temp_dic = record.pos, self.forward_dic
            corr_depth = self._get_gc_corr_dep(record.pos)
        try:
            temp_dic[curr_pos] += corr_depth
        except KeyError:
            temp_dic[curr_pos] = corr_depth
            self._call_output(temp_dic)

    def writetofile(self):
        with open(self.arg.out, 'w') as f_output:
            for key, value in self.outputdic.iteritems():
                f_output.write('{}\t{}\n'.format(key, value))

    def call(self):
        self._call_output(self.forward_dic, max_lst_range=0, max_size=0)
        self._call_output(self.reverse_dic, max_lst_range=0, max_size=0)

    def reset(self):
        self.forward_dic = {}
        self.reverse_dic = {}

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
        fasta_str = self._fasta.fetch_string(self.chrom, pos,
                                             self._GC_model_len)
        gc_idx = fasta_str.count('G')+fasta_str.count('C')
        return self._model[gc_idx]


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('--gcmodel', help="fastafile")
    parser.add_argument('--bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_phasogram.txt')
    return parser.parse_args(argv)


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')[:3]
                chrom = input_line.pop(0).replace('chr', '')
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        yield (args.chrom, args.start, args.end)


def main(argv):
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    phaso = Phasogram(args)
    for chrom, start, end in read_bed(args):
        last_tid = None
        phaso.reset()
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            if last_tid and record.tid != last_tid:
                phaso.call()
                phaso.reset()
            phaso.update(record, chrom)
            last_tid = record.tid
        phaso.call()
    phaso.writetofile()
    samfile.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
