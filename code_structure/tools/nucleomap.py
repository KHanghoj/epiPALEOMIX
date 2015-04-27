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
    read_bed_WO, strtobool, GC_correction, corr_fasta_chr  # , Cache
import gzip


class Nucleosome_Prediction(GC_correction):
    """docstring for Nucleosome_Prediction"""
    def __init__(self, arg):
        self.arg = arg
        self._SIZE = self.arg.SIZE
        self._OFFSET, self._NEIGHBOR = self.arg.OFFSET, self.arg.FLANKS
        self._POSITION_OFFSET = self._OFFSET+self._NEIGHBOR
        self._TOTAL_WIN_LENGTH = self._SIZE+(2*self._OFFSET)+(2*self._NEIGHBOR)
        self._CENTERINDEX = (self._SIZE-1)/2
        GC_correction.__init__(self)
        self._mindepth = int(self.arg.MinDepth)
        self._seq_len = int(self.arg.DequeLen)
        self._zeros = [0]*self._seq_len
        self._read_max_len = self._seq_len-(self.arg.MaxReadLen+100)
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
                nwisedat = self._nwise_cut()
                self._call_window(nwisedat)  # call up to actual idx
                self._deq_depth.extend([0]*(self._actual_idx -
                                            self._TOTAL_WIN_LENGTH))
            self._actual_idx = self._TOTAL_WIN_LENGTH
            self._last_ini = record.pos-self._TOTAL_WIN_LENGTH
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

    def _nwise_cut(self):
        return self._nwise_zip(islice(self._deq_depth, 0, self._actual_idx-1),
                               islice(self._deq_pos, 0, self._actual_idx-1))

    def _nwise(self, deque_lst, n=None):
        n = self._TOTAL_WIN_LENGTH if n is None else n
        return izip(*(islice(g, i, None)
                      for i, g in enumerate(tee(deque_lst, n))))

    def _nwise_zip(self, deq, pos):
        return izip(self._nwise(deq), self._nwise(pos))

    def call_final_window(self):
        self._call_window(self._nwise_zip(self._deq_depth, self._deq_pos))
        
    def _call_window(self, nwisedat):
        ''' docstring '''
        for self.win_depth, self.win_position in nwisedat:
            # returns two tuples that are unpacked by
            # self.win_depth and self.win_position
            if (self.win_depth[self._POSITION_OFFSET+self._CENTERINDEX] <
                self._mindepth):
                continue
            dep_posses = self._call_max(self.win_depth[self._POSITION_OFFSET:
                                                       self._POSITION_OFFSET +
                                                       self._SIZE])
            if dep_posses:
                spacerL = (sum(self.win_depth[:self._NEIGHBOR]) /
                           float(self._NEIGHBOR))
                spacerR = (sum(self.win_depth[-self._NEIGHBOR:]) /
                           float(self._NEIGHBOR))
                # both cannot be 0
                if spacerL or spacerR:
                    center_depth, min_idx, max_idx = dep_posses
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
            val = self.window[(self._CENTERINDEX + idx)]
            if val != self.maxdepth:
                break
            self.call.append(self._CENTERINDEX+idx)

    def _call_max(self, window):
        ''' docstring '''
        self.call = list()
        self.window = window
        self.maxdepth = max(self.window)
        if self.maxdepth >= self._mindepth:
            if self.window[self._CENTERINDEX] == self.maxdepth:
                self.call.append(self._CENTERINDEX)
                self._check_width(1, self._CENTERINDEX, 1)
                self._check_width(-1, -self._CENTERINDEX, -1)
                return (self.maxdepth, self._POSITION_OFFSET+min(self.call),
                        self._POSITION_OFFSET+max(self.call))

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

    def closefile(self):
        self.f_output.close()
        self._fasta_dat.closefile()

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
    parser.add_argument('--NucleosomeSize', dest='SIZE', help="..", type=int, default=147)
    parser.add_argument('--NucleosomeFlanks', dest='FLANKS', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeOffset', dest='OFFSET', help="..", type=int, default=12)
    return parser.parse_known_args(argv)


def run(args):
    read_bed = read_bed_W if strtobool(args.BamChromType) else read_bed_WO
    args.FastaChromType = strtobool(args.FastaChromType) if args.GCmodel else True
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
