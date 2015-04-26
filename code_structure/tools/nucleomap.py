#!/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
from collections import deque
from itertools import islice, izip, tee
from os.path import exists, splitext
from shutil import move
from epiomix_commonutils import read_bed_W, \
    read_bed_WO, strtobool, Cache, GC_correction, corr_fasta_chr
import gzip
# CONSTANTS:
_SIZE = 147  # the window
_OFFSET = 12  # _OFFSET between spacer and nucleosomal DNA
_NEIGHBOR = 25  # flanking regions to be considered for log-odd ration score
_POSITION_OFFSET = _NEIGHBOR+_OFFSET
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_CENTERINDEX = (_SIZE-1)/2


class Nucleosome_Prediction(GC_correction):
    """docstring for Nucleosome_Prediction"""
    def __init__(self, arg):
        self.arg = arg
        self._fasta = Cache(self.arg.FastaPath)
        GC_correction.__init__(self)
        self._mindepth = int(self.arg.MinDepth)
        self._seq_len = int(self.arg.DequeLen)
        self._zeros = [0]*self._seq_len
        self._read_max_len = self._seq_len-(self.arg.MaxReadLen+50)
        self._deq_depth = deque(maxlen=self._seq_len)
        self._deq_pos = deque(maxlen=self._seq_len)
        self._last_ini = -(self._seq_len+1)
        self._actual_idx, self.f_output = None, None
        self._GC_model_len = 0
        self._output_dic = dict()
        self._fmt = '{0}\t{1}\t{2}\t{3}\t{4}\t{bedcoord}\n'
        self._GCmodel_ini()
        self._makeoutputfile()

    def update_depth(self, record):
        self._actual_idx = (record.pos - self._last_ini)
        if self._actual_idx > self._read_max_len:
            if self._actual_idx > self._seq_len:
                nwisedat = self._nwise_zip(self._deq_depth, self._deq_pos)
                self._call_window(nwisedat)  # call everything
                self._deq_depth.extend(self._zeros)
                self.writetofile()
            else:
                tmp_deq_dep = islice(self._deq_depth, 0, self._actual_idx-1)
                tmp_deq_pos = islice(self._deq_pos, 0, self._actual_idx-1)
                nwisedat = self._nwise_zip(tmp_deq_dep,
                                           tmp_deq_pos)
                self._call_window(nwisedat)  # call up to actual idx
                self._deq_depth.extend([0]*(self._actual_idx -
                                            _TOTAL_WIN_LENGTH))
            self._actual_idx = _TOTAL_WIN_LENGTH
            self._last_ini = record.pos-_TOTAL_WIN_LENGTH
            self._deq_pos.extend((pos+1) for
                                 pos in xrange(self._last_ini, self._last_ini +
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

    def _nwise(self, deque_lst, n=_TOTAL_WIN_LENGTH):
        return izip(*(islice(g, i, None)
                      for i, g in enumerate(tee(deque_lst, n))))

    def _nwise_zip(self, deq, pos):
        return izip(self._nwise(deq), self._nwise(pos))

    def call_final_window(self):
        nwisedat = self._nwise_zip(self._deq_depth, self._deq_pos)
        self._call_window(nwisedat)

    def _call_window(self, nwisedat):
        ''' docstring '''
        for self.win_depth, self.win_position in nwisedat:
            # returns two tuples that are unpacked by
            # self.win_depth and self.win_position
            if self.win_depth[_POSITION_OFFSET+_CENTERINDEX] < self._mindepth:
                continue
            dep_min_max_pos = self._call_max(self.win_depth[_POSITION_OFFSET:
                                                            _POSITION_OFFSET +
                                                            _SIZE])
            if dep_min_max_pos:
                spacerL = sum(self.win_depth[:_NEIGHBOR])/float(_NEIGHBOR)
                spacerR = sum(self.win_depth[-_NEIGHBOR:])/float(_NEIGHBOR)
                # both cannot be 0
                if spacerL or spacerR:
                    center_depth, min_idx, max_idx = dep_min_max_pos
                    sizeofwindow = (max_idx-min_idx)
                    mean_spacer = 1.0 + (0.5 * (spacerL + spacerR))
                    # dividing by width of nucleosome called
                    score = ((float(center_depth) / mean_spacer) /
                             (sizeofwindow+1.0))

                    start_pos = self.win_position[min_idx]
                    end_pos = self.win_position[max_idx]
                    self._output_dic[start_pos] = [end_pos, center_depth, score]

    def _check_width(self, start, end, incre):
        for idx in xrange(start, end, incre):
            val = self.window[(_CENTERINDEX + idx)]
            if val != self.maxdepth:
                break
            self.call.append(_CENTERINDEX+idx)

    def _call_max(self, window):
        ''' docstring '''
        self.call = list()
        self.window = window
        self.maxdepth = max(self.window)
        if self.maxdepth >= self._mindepth:
            if self.window[_CENTERINDEX] == self.maxdepth:
                self.call.append(_CENTERINDEX)
                self._check_width(1, _CENTERINDEX, 1)
                self._check_width(-1, -_CENTERINDEX, -1)
                return (self.maxdepth, _POSITION_OFFSET+min(self.call),
                        _POSITION_OFFSET+max(self.call))

    def writetofile(self):
        ''' dfs '''
        if self._output_dic:
            for key, val in iter(sorted(self._output_dic.iteritems())):
                self.f_output.write(self._fmt.format(self.chrom,
                                    key, *val, bedcoord=self.bedcoord))
            self._output_dic.clear()

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
        # header = '#chrom\tstart\tend\tdepth\tscore\tbedcoord\n'
        # self.f_output.write(header)

    def closefile(self):
        self.f_output.close()
        self._fasta.closefile()

    def reset_deques(self, chrom, start, end):
        self._deq_depth.clear()
        self._deq_pos.clear()
        self._output_dic.clear()
        self._last_ini = -(self._seq_len+1)
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = '{}_{}_{}'.format(self.chrom, self.start, self.end)


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--MaxReadLen', help="..", type=int, default=150)
    parser.add_argument('--DequeLen', help="..", type=int, default=2000)
    parser.add_argument('--MinDepth', help="..", type=int, default=5)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--GCmodel', help='..', type=str, default=None)
    parser.add_argument('--FastaChromType', dest='FastaChromType')
    parser.add_argument('--BamChromType', dest='BamChromType')
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    return parser.parse_known_args(argv)


def run(args):
    read_bed = read_bed_W if strtobool(args.BamChromType) else read_bed_WO
    args.FastaChromType = strtobool(args.FastaChromType)
    samfile = pysam.Samfile(args.bam, "rb")
    nucl_pred_cls = Nucleosome_Prediction(args)
    for chrom, start, end in read_bed(args):
        nucl_pred_cls.reset_deques(corr_fasta_chr(args, chrom), start, end)
        for record in samfile.fetch(chrom, start, end):
            nucl_pred_cls.update_depth(record)
        nucl_pred_cls.call_final_window()
        nucl_pred_cls.writetofile()
    nucl_pred_cls.closefile()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
