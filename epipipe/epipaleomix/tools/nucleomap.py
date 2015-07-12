#!/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
from collections import deque, namedtuple
from os.path import exists, splitext
from shutil import move
from epipaleomix.tools.commonutils import \
    GC_correction, read_bed
import gzip


class Nucleosome_Prediction(GC_correction):
    """docstring for Nucleosome_Prediction"""
    def __init__(self, arg):
        self.arg = arg
        self._SIZE = self.arg.SIZE
        self._OFFSET, self._NEIGHBOR = self.arg.OFFSET, self.arg.FLANKS
        self._POSITION_OFFSET = self._OFFSET+self._NEIGHBOR
        self._TOTAL_WIN_LENGTH = self._SIZE+(2*self._OFFSET)+(2*self._NEIGHBOR)
        self._HALFWINDOW =  (self._TOTAL_WIN_LENGTH/2)+1
        self._CENTERINDEX = (self._SIZE-1)/2
        GC_correction.__init__(self)
        self._mindepth = int(self.arg.MinDepth)
        self._seq_len = int(self._TOTAL_WIN_LENGTH*2)
        self._zeros = [0]*self._seq_len
        self._deq_depth = deque(self._zeros, maxlen=self._seq_len)
        self._mainlist = []
        self._last_ini = None
        self.f_output = None
        self._GC_model_len = 0
        self._outputlist = []
        self._fmt = '{}\t{}\t{}\t{}\t{}\t{bedcoord}\n'
        self._GCmodel_ini()
        self._makeoutputfile()
        self._calltuple = namedtuple('call', 's e depth score')

    def update_depth(self, record):
        if not self._last_ini:
            self._last_ini = record.pos+1
            self._last_pos = record.pos
        _jump_idx = record.pos - self._last_pos
        if _jump_idx:
            if _jump_idx >= self._seq_len:
                self._mainlist.extend(self._deq_depth)
                self._call_window()
                self._mainlist = self._zeros[:self._HALFWINDOW]
                self._last_ini = record.pos+1-self._HALFWINDOW
                self._deq_depth = deque(self._zeros, maxlen=self._seq_len)
                _jump_idx = 0
            else:
                while _jump_idx:
                    self._mainlist.append(self._deq_depth.popleft())
                    self._deq_depth.append(0)
                    _jump_idx -= 1
        self._last_pos = record.pos

        if record.is_reverse:
            corr_depth = self._get_gc_corr_dep(record.aend-self._GC_model_len)
        else:
            corr_depth = self._get_gc_corr_dep(record.pos)

        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(_jump_idx, _jump_idx + count):
                    self._deq_depth[idx] += corr_depth
                _jump_idx += count
            elif cigar in (2, 3, 6):
                _jump_idx += count

    def _depthwindows(self):
        n = len(self._mainlist)
        for idx in xrange(0, n-self._TOTAL_WIN_LENGTH+1):
            yield idx, self._mainlist[idx:(idx+self._TOTAL_WIN_LENGTH)]
    
    def call_final_window(self):
        self._mainlist.extend(self._deq_depth)
        self._call_window()
        
    def _call_window(self):
        ''' docstring '''
        lasttup = ()
        for idx, win_depth in self._depthwindows():
            if win_depth[self._POSITION_OFFSET+self._CENTERINDEX] < self._mindepth:
                continue
            dep_posses = self._call_max(win_depth[self._POSITION_OFFSET:
                                                  self._POSITION_OFFSET +
                                                  self._SIZE])
            if dep_posses:
                spacerL = (sum(win_depth[:self._NEIGHBOR]) /
                           float(self._NEIGHBOR))
                spacerR = (sum(win_depth[-self._NEIGHBOR:]) /
                           float(self._NEIGHBOR))
                # both cannot be 0
                if spacerL or spacerR:
                    center_depth, min_idx, max_idx = dep_posses
                    sizeofwindow = (max_idx-min_idx)
                    mean_spacer = (0.5 * (spacerL + spacerR))
                    if mean_spacer < 1: # when GCcorrection, we see very low flanks. gives unrealistic high score
                        mean_spacer = 1
                    # divide peak by mean of flanks then by width of nucleosome
                    score = ((float(center_depth) / mean_spacer) /
                             (sizeofwindow+1.0))
                    start_pos = idx+self._last_ini+min_idx
                    end_pos = idx+self._last_ini+max_idx

                    if start_pos >= self.start and start_pos <= self.end:
                        if not lasttup: # initialize
                            lasttup = self._calltuple(start_pos, end_pos, center_depth, score)
                            continue
                        if start_pos <= lasttup.e: # the nucl dyad/center overlap
                            if score > lasttup.score:
                                lasttup = self._calltuple(start_pos, end_pos, center_depth, score)
                        else:  # they do not overlap. send lasttup to outputlist
                            self._outputlist.append(lasttup)
                            lasttup = self._calltuple(start_pos, end_pos, center_depth, score)
        if lasttup:
            self._outputlist.append(lasttup)

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
        if self._outputlist:
            for dat in self._outputlist: 
                self.f_output.write(self._fmt.format(self.chrom,
                                                     *dat,
                                                     bedcoord=self.bedcoord))
            self._outputlist = []

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')

    def closefile(self):
        self.f_output.close()
        try:
            self._fasta_dat.closefile()
        except AttributeError:
            pass

    def reset_deques(self, chrom, start, end, bedcoord):
        self._deq_depth = deque(self._zeros, maxlen=self._seq_len)
        self._outputlist = []
        self._last_ini = None
        self._mainlist = [] # everytime new bed is called. flanks are made automatically. see run func
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = bedcoord


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--MinDepth', help="..", type=int, default=5)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--GCmodel', help='..', type=str, default=None)
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeSize', dest='SIZE', help="..", type=int, default=147)
    parser.add_argument('--NucleosomeFlanks', dest='FLANKS', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeOffset', dest='OFFSET', help="..", type=int, default=12)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.Samfile(args.bam, "rb")
    nucl_pred_cls = Nucleosome_Prediction(args)
    flanks = (nucl_pred_cls._TOTAL_WIN_LENGTH/2)+1
    for chrom, start, end, bedcoord in read_bed(args):
        nucl_pred_cls.reset_deques(chrom, start, end, bedcoord)
        start = 0 if start-flanks < 0 else start-flanks
        for record in samfile.fetch(chrom, start, end+flanks):
            if record.mapq < args.MinMappingQuality:
                continue  # do not analyze low quality records
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
