#!/opt/local/bin/python
from __future__ import print_function

#  from fileinput import inpua
import sys
import pysam
import math
import argparse
from collections import deque
from itertools import islice, izip, tee
from os import remove
#### Constants:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_POSITION_OFFSET = _NEIGHBOR+_OFFSET
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_CENTERINDEX = (_SIZE-1)/2
# _MAXLEN = int(4e2)+110
# _MAXLEN = int(1e3)+110
_MAXLEN = 500+110
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

    def fetch_string(self, chrom, start, nbases=2):
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
    def __init__(self, samfile, outputpath, seq_len=_MAXLEN):
        # if new chrom: reset every thing in main():
        self.samf = samfile
        self.outputpath = outputpath
        self._seq_len = int(seq_len)
        self._deq_depth = deque(maxlen=self._seq_len)
        self._deq_pos = deque(maxlen=self._seq_len)
        self._last_ini = -self._seq_len
        self._actual_idx = None
        self._present_tid = ''
        self._zeros = [0]*self._seq_len
        self._output_dic = dict()

    def update_depth(self, record):
        self.jump = (record.pos-self._last_ini)-self._seq_len
        # if self.jump+record.alen >= 0:
        if self.jump+record.alen >= -110:
            # it could be that we need to let the window slide to the end
            # that means we need to add 110 zeros at the end
            self.call_window()

            if self.jump > _TOTAL_WIN_LENGTH:
                # next read is further away than entire window
                self._deq_depth.extend(self._zeros)
                if len(self._output_dic) > 1000:
                    self.writetofile()
            else:
                self._deq_depth.extend(self._zeros[:self._seq_len -
                                                   _TOTAL_WIN_LENGTH +
                                                   self.jump])
            self._last_ini = record.pos-_TOTAL_WIN_LENGTH
            self._deq_pos.extend((pos+1) for pos in xrange(self._last_ini,
                                                           self._last_ini +
                                                           self._seq_len))
        self._present_tid = record.tid
        self._actual_idx = (record.pos-self._last_ini)
        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(self._actual_idx, self._actual_idx + count):
                    self._deq_depth[idx] += 1
                self._actual_idx += count
            elif cigar in (2, 3, 6):
                self._actual_idx += count

    def move_deqs(self, record, var1):
        # this is not being used
        self._last_ini = record.pos-((_TOTAL_WIN_LENGTH-1)/var1)
        pass

    def _nwise(self, n=_TOTAL_WIN_LENGTH):
        self.depths = izip(*(islice(g, i, None)
                           for i, g in enumerate(tee(self._deq_depth, n))))
        self.posis = izip(*(islice(g, i, None)
                            for i, g in enumerate(tee(self._deq_pos, n))))
        return izip(self.depths, self.posis)

    def call_window(self):
        ''' docstring '''
        for self.win_depth, self.win_positions in self._nwise():
            # if islice(self.win_depth,
            #           _POSITION_OFFSET+_CENTERINDEX, None).next() <= _MINDEPTH:
            #     continue
            if self.win_depth[_POSITION_OFFSET+_CENTERINDEX] <= _MINDEPTH:
                continue
            calls_center = self._call_max(self.win_depth[_POSITION_OFFSET:
                                          _POSITION_OFFSET+_SIZE])
            # calls_center = self._call_max(islice(self.win_depth,
            #                               _POSITION_OFFSET,
            #                               _POSITION_OFFSET+_SIZE))

            if calls_center:

                # spacerL = islice(self.win_depth, 0, _NEIGHBOR)
                # spacerR = islice(self.win_depth,
                #                  _POSITION_OFFSET+_SIZE+_OFFSET, None)
                spacerL = self.win_depth[:_NEIGHBOR]
                spacerR = self.win_depth[-_NEIGHBOR:]
                calls_spacerL = sum(spacerL)/float(_NEIGHBOR)
                calls_spacerR = sum(spacerR)/float(_NEIGHBOR)
                mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)
                center_depth = calls_center.itervalues().next()
                score = float(center_depth) / mean_spacer
                min_pos = _POSITION_OFFSET+min(calls_center)
                max_pos = _POSITION_OFFSET+max(calls_center)
                start_pos = self.win_positions[min_pos]
                end_pos = self.win_positions[max_pos]
                # start_pos = islice(self.win_positions, min_pos, None).next()
                # end_pos = islice(self.win_positions, max_pos, None).next()

                key = start_pos
                if key in self._output_dic.keys():
                    if self._output_dic[key][-1] < score:
                        self._output_dic[key] = [end_pos, center_depth, score]
                else:
                    self._output_dic[key] = [end_pos, center_depth, score]

    def _check_width(self, start, end, incre):
        for idx in xrange(start, end, incre):
            val = self.window[(_CENTERINDEX + idx)]
            if val != self.maxdepth:
                break
            self.call[(_CENTERINDEX+idx)] = val

    def _call_max(self, lst):
        ''' docstring '''
        self.call = {}  # this is the return dic with indexposition and depth
        self.window = lst
        self.maxdepth = max(self.window)
        if self.maxdepth > _MINDEPTH:
            if self.window[_CENTERINDEX] == self.maxdepth:
                self.call[_CENTERINDEX] = self.window[_CENTERINDEX]
                self._check_width(1, _CENTERINDEX, 1)
                self._check_width(-1, -_CENTERINDEX, -1)
                return self.call

    def writetofile(self):
        ''' dfs '''
        # append to file instead
        if self._output_dic:
            chrom = self.samf.getrname(self._present_tid)
            key_values = iter(sorted([key]+value for key,
                              value in
                              self._output_dic.iteritems()))
            fmt = '{0}\t{1}\t{2}\t{3}\t{4}\n'
            for dat in key_values:
                self.f_output.write(fmt.format(chrom, *dat))
            self._output_dic.clear()

    def makeoutputfile(self):
        # i want to write to file every chrom, to speed up:
        try:
            remove(self.outputpath)
            self.f_output = open(self.outputpath, 'a')
        except OSError:
            self.f_output = open(self.outputpath, 'a')

    def closefile(self):
        return self.f_output.close()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--bed', help="...")
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
    nucl_pred_cls = nucleosome_prediction(samfile, args.out)
    nucl_pred_cls.makeoutputfile()
    for chrom, start, end in read_bed(args, ''):
        last_tid = ''
        for record in samfile.fetch(chrom, start, end):
            if record.tid != last_tid:
                nucl_pred_cls.writetofile()
            last_tid = record.tid
            nucl_pred_cls.update_depth(record)
        nucl_pred_cls.call_window()
        nucl_pred_cls.writetofile()
    nucl_pred_cls.writetofile()
    nucl_pred_cls.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
