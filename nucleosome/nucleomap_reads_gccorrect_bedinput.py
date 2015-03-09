#!/opt/local/bin/python
from __future__ import print_function

#  from fileinput import inpua
import sys
import pysam
import argparse
from itertools import islice, izip, tee
# CONSTANTS:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_POSITION_OFFSET = _NEIGHBOR+_OFFSET
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_CENTERINDEX = (_SIZE-1)/2
_HALF_WIN_LENGTH = (_TOTAL_WIN_LENGTH-1)/2
_MAXLEN = 500
_READ_MAX_LEN = _MAXLEN-_HALF_WIN_LENGTH
_MINDEPTH = 20


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


class nucleosome_prediction(object):
    """docstring for nucleosome_prediction"""
    def __init__(self, arg):
        # if new chrom: reset every thing in main():
        self.arg = arg
        self.fasta = Cache(self.arg.fastafile)
        self._outputpath = self.arg.out
        self._deq_depth = list()
        self._deq_pos = list()
        self._last_ini = None
        self._actual_idx = None
        self._output_dic = dict()

    def update_depth(self, record, chrom):
        self._present_chrom = chrom
        self._actual_idx = (record.pos-self._last_ini)+_TOTAL_WIN_LENGTH
        if record.is_reverse:
            # gc_idx = self.get_gc_count(record.aend-len(self.model))
            gc_idx = self.get_gc_count(record.aend-self._GC_model_len)
        else:
            gc_idx = self.get_gc_count(record.pos+1)
        corr_depth = (self.model[gc_idx])
        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(self._actual_idx, self._actual_idx + count):
                    self._deq_depth[idx] += corr_depth
                self._actual_idx += count
            elif cigar in (2, 3, 6):
                self._actual_idx += count

    def _nwise(self, n=_TOTAL_WIN_LENGTH):
        self.depths = izip(*(islice(g, i, None)
                           for i, g in enumerate(tee(self._deq_depth, n))))
        self.posis = izip(*(islice(g, i, None)
                            for i, g in enumerate(tee(self._deq_pos, n))))
        # izip returns two tuples for each entry zipped together.
        # unpacked by two variables in for loop in call_window
        return izip(self.depths, self.posis)

    def call_window(self):
        ''' docstring '''
        for self.win_depth, self.win_positions in self._nwise():
            # returns two tuples that are unpacked by
            # self.win_depth and self.win_positions
            if self.win_depth[_POSITION_OFFSET+_CENTERINDEX] <= _MINDEPTH:
                continue
            calls_center = self._call_max(self.win_depth[_POSITION_OFFSET:
                                          _POSITION_OFFSET+_SIZE])
            # here calls center is just a number. the max depth
            if calls_center:
                spacerL = self.win_depth[:_NEIGHBOR]
                spacerR = self.win_depth[-_NEIGHBOR:]
                calls_spacerL = sum(spacerL)/float(_NEIGHBOR)
                calls_spacerR = sum(spacerR)/float(_NEIGHBOR)
                mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)
                score = float(calls_center) / mean_spacer
                key = self.win_positions[_POSITION_OFFSET+_CENTERINDEX]
                if score > 1:
                    try:
                        self._output_dic[key] += 1
                    except KeyError:
                        self._output_dic[key] = 1

    def _call_max(self, lst):
        ''' docstring '''
        self.window = lst
        self.maxdepth = max(self.window)
        if self.maxdepth > _MINDEPTH:
            if self.window[_CENTERINDEX] == self.maxdepth:
                if self.window[_CENTERINDEX-1] == self.maxdepth or \
                        self.window[_CENTERINDEX+1] == self.maxdepth:
                    return None
                return self.maxdepth

    def writetofile(self):
        ''' dfs '''
        # append to file instead
        with open(self._outputpath, 'w') as f:
            mininum = min(self._output_dic.iterkeys())
            maximum = max(self._output_dic.iterkeys())
            fmt = '{0}\t{1}\n'
            for key in xrange(mininum, maximum+1, 1):
                value = self._output_dic.get(key, 0)
                f.write(fmt.format(key, value))

    def reset_deques(self, start, end):
        deq_length = end-start
        self._deq_depth = list([0.0]*(deq_length+(2*_TOTAL_WIN_LENGTH)))
        self._deq_pos = list(xrange(-_TOTAL_WIN_LENGTH,
                                    deq_length+_TOTAL_WIN_LENGTH))
        self._last_ini = start

    def gcmodel_ini(self):
        if self.arg.gcmodel:
            with open(self.arg.gcmodel, 'r') as f:
                self.model = [float(line.rstrip('\n').split('\t')[-1])
                              for line in f]
                self._GC_model_len = len(self.model)
        else:
            self._GC_model_len = 0  # this is default if not assign
            self.model = [1]*(self._GC_model_len+1)

    def get_gc_count(self, pos):
        fasta_str = self.fasta.fetch_string(self._present_chrom, pos,
                                            self._GC_model_len)
        return(fasta_str.count('G')+fasta_str.count('C'))


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="...")
    parser.add_argument('--gcmodel', help="pathtogcfile", default=None)
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", default=None)
    parser.add_argument('--end', help="...", default=None)
    parser.add_argument('--out', help='...', default='out_nucleomap.txt')
    return parser.parse_args(argv)


def read_bed(args, chromtype):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = line.rstrip('\n').split('\t')
                chrom = input_line.pop(0).replace('chr', chromtype)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        try:
            yield (args.chrom, int(args.start), int(args.end))
        except TypeError:
            yield (args.chrom, args.start, args.end)


def main(argv):
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    nucl_pred_cls = nucleosome_prediction(args)
    nucl_pred_cls.gcmodel()
    for chrom, start, end in read_bed(args, ''):
        nucl_pred_cls.reset_deques(start, end)
        for record in samfile.fetch(chrom, start, end):
            nucl_pred_cls.update_depth(record, chrom)
        nucl_pred_cls.call_window()
    nucl_pred_cls.writetofile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
