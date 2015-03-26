#!/opt/local/bin/python
from __future__ import print_function
# python ~/research/projects/epiomix/nucleosome/nucleomap_reads_gccorrect.py
# chrome.fa test.bam --chrom 22 --start 16000000 --end 17000000 --out hmm.txt

# python ~/research/projects/epiomix/nucleosome/nucleomap_reads_gccorrect.py
# chrome.fa test.bam --gcmodel gccorrect_saq.txt --chrom 22 --start 16000000
# --end 17000000 --out hmm.txt
import sys
import pysam
import argparse
from collections import deque
from itertools import islice, izip, tee
# CONSTANTS:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_POSITION_OFFSET = _NEIGHBOR+_OFFSET
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_CENTERINDEX = (_SIZE-1)/2
# _HALF_WIN_LENGTH = (_TOTAL_WIN_LENGTH-1)/2
_MAXLEN = 1000
# _READ_MAX_LEN = _MAXLEN-(_TOTAL_WIN_LENGTH*2)
_MINDEPTH = 4


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


class Nucleosome_Prediction(object):
    """docstring for Nucleosome_Prediction"""
    def __init__(self, arg, seq_len=_MAXLEN, mindepth=_MINDEPTH):
        self.arg = arg
        self._fasta = Cache(self.arg.fastafile)
        self._mindepth, self._seq_len = int(mindepth), int(seq_len)
        self._zeros = [0]*(self._seq_len+1)
        self._read_max_len = self._seq_len-(_TOTAL_WIN_LENGTH)
        self._deq_depth = deque(maxlen=self._seq_len)
        self._deq_pos = deque(maxlen=self._seq_len)
        self._last_ini = -self._seq_len
        self._actual_idx = None
        self._output_dic = dict()
        self._gcmodel_ini()
        # self._makeoutputfile()

    def update_depth(self, record):
        self._actual_idx = (record.pos - self._last_ini)
        if self._actual_idx >= self._read_max_len:
            if self._actual_idx > self._seq_len+_TOTAL_WIN_LENGTH:
                self._nwise(self._deq_depth, self._deq_pos)
                self._call_window()  # call everything
                self._deq_depth.extend(self._zeros)
                # self.writetofile()
            else:
                self._nwise(islice(self._deq_depth,
                                   0, self._actual_idx-1),  # MINUS 1 FIX EVERY
                            islice(self._deq_pos,
                                   0, self._actual_idx-1))
                self._call_window()  # call up to actual idx
                self._deq_depth.extend(self._zeros[:self._actual_idx -
                                                   _TOTAL_WIN_LENGTH])

            self._actual_idx = _TOTAL_WIN_LENGTH
            self._last_ini = record.pos-_TOTAL_WIN_LENGTH
            self._deq_pos.extend((pos+1) for pos in xrange(self._last_ini,
                                                           self._last_ini +
                                                           self._seq_len))
        if record.is_reverse:
            corr_depth = self._get_gc_corr_dep(record.aend-self._GC_model_len)
        else:
            corr_depth = self._get_gc_corr_dep(record.pos)
        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(self._actual_idx, self._actual_idx + count):
                    self._deq_depth[idx] += corr_depth
                self._actual_idx += count
            elif cigar in (2, 3, 6):
                self._actual_idx += count

    def _nwise(self, dep, pos, n=_TOTAL_WIN_LENGTH):
        '''izip returns two tuples for each entry zipped together.
        unpacked by two variables in for loop in _call_window '''
        depths = izip(*(islice(g, i, None)
                      for i, g in enumerate(tee(dep, n))))
        posis = izip(*(islice(g, i, None)
                     for i, g in enumerate(tee(pos, n))))
        self.nwise_dat = izip(depths, posis)

    # # this is new and not called:
    # ######
    # def _nwise1(self, deque_lst, n=_TOTAL_WIN_LENGTH):
    #     '''izip returns two tuples for each entry zipped together.
    #     unpacked by two variables in for loop in _call_window '''
    #     return izip(*(islice(g, i, None)
    #                   for i, g in enumerate(tee(deque_lst, n))))

    # def nwise_zip(deq, pos):
    #     self.nwise_dat = izip(self._nwise1(deq), self._nwise1(pos))
    # ######

    def call_final_window(self):
        self._nwise(self._deq_depth, self._deq_pos)
        self._call_window()

    def _call_window(self):
        ''' docstring '''
        for self.win_depth, self.win_position in self.nwise_dat:
            # returns two tuples that are unpacked by
            # self.win_depth and self.win_position
            if self.win_depth[_POSITION_OFFSET+_CENTERINDEX] <= self._mindepth:
                continue
            dep_min_max_pos = self._call_max(self.win_depth[_POSITION_OFFSET:
                                                            _POSITION_OFFSET +
                                                            _SIZE])
            if dep_min_max_pos:
                center_depth, min_idx, max_idx = dep_min_max_pos
                if min_idx != max_idx:
                    continue
                spacerL = sum(self.win_depth[:_NEIGHBOR])/float(_NEIGHBOR)
                spacerR = sum(self.win_depth[-_NEIGHBOR:])/float(_NEIGHBOR)
                mean_spacer = 1.0 + 0.5 * (spacerL + spacerR)
                score = float(center_depth) / mean_spacer
                centerofwin = (max_idx-min_idx)
                score = score/(centerofwin+1.0)
                key = self.win_position[(min_idx+((centerofwin)/2))]
                if score > 1:
                    try:
                        self._output_dic[key-self.start-1] += 1
                    except KeyError:
                        self._output_dic[key-self.start-1] = 1

    def _check_width(self, start, end, incre):
        for idx in xrange(start, end, incre):
            val = self.window[(_CENTERINDEX + idx)]
            if val != self.maxdepth:
                break
            self.call.append(_CENTERINDEX+idx)

    def _call_max(self, window):
        ''' docstring '''
        self.call = list()  # this is the return dic with indexposition and depth
        self.window = window
        self.maxdepth = max(self.window)
        if self.maxdepth > self._mindepth:
            if self.window[_CENTERINDEX] == self.maxdepth:
                self.call.append(_CENTERINDEX)
                self._check_width(1, _CENTERINDEX, 1)
                self._check_width(-1, -_CENTERINDEX, -1)
                return (self.maxdepth, _POSITION_OFFSET+min(self.call),
                        _POSITION_OFFSET+max(self.call))

    # def writetofile(self):
    #     ''' dfs '''
    #     if self._output_dic:
    #         fmt = '{0}\t{1}\t{2}\t{3}\t{4}\n'
    #         for key, val in iter(sorted(self._output_dic.iteritems())):
    #             self.f_output.write(fmt.format(self.chrom, key, *val))
    #         self._output_dic.clear()
    #         # NOTE OUTPUTDIC IS CLEARED/EMPTIED
    #         # AFTER WRITETOFILE

    def writetofile(self):
        ''' dfs '''
        # append to file instead
        with open(self.arg.out, 'w') as f:
            mininum = min(self._output_dic.iterkeys())
            maximum = max(self._output_dic.iterkeys())
            fmt = '{0}\t{1}\n'
            for key in xrange(mininum, maximum+1, 1):
                value = self._output_dic.get(key, 0)
                f.write(fmt.format(key, value))
        self._fasta.closefile()

    # def _makeoutputfile(self):
    #     ''' want to write to file every chrom, to keep scalablility'''
    #     try:
    #         remove(self._outputpath)
    #         self.f_output = open(self._outputpath, 'a')
    #     except OSError:
    #         self.f_output = open(self._outputpath, 'a')

    # def closefile(self):
    #     return self.f_output.close()

    def reset_deques(self, chrom, start, end):
        self._deq_depth.clear()
        self._deq_pos.clear()
        self._last_ini = -self._seq_len
        self.chrom, self.start, self.end = chrom, start, end

    def _gcmodel_ini(self):
        self._model = list()
        if self.arg.gcmodel:
            with open(self.arg.gcmodel, 'r') as f:
                self._model = [float(line.rstrip('\n').split('\t')[-1])
                               for line in f]
                self._GC_model_len = len(self._model)
        else:
            self._GC_model_len = 0  # this is default if not assign
            self._model = [1]*(self._GC_model_len+1)

    def _get_gc_corr_dep(self, pos):
        if self._GC_model_len:
            fasta_str = self._fasta.fetch_string(self.chrom,
                                                 pos, self._GC_model_len)
            gc_idx = fasta_str.count('G')+fasta_str.count('C')
            return self._model[gc_idx]
        else:
            return 1


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="...")
    parser.add_argument('--gcmodel', help="fastafile")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", default=None)
    parser.add_argument('--end', help="...", default=None)
    parser.add_argument('--out', help='...', default='out_nucleomap.txt')
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
    nucl_pred_cls = Nucleosome_Prediction(args)
    for chrom, start, end in read_bed(args):
        # last_tid = ''
        nucl_pred_cls.reset_deques(chrom, start, end)
        for record in samfile.fetch(chrom, start, end):
            nucl_pred_cls.update_depth(record)
        nucl_pred_cls.call_final_window()
    nucl_pred_cls.writetofile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
