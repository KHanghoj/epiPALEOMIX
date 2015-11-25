#!/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
import math
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
        self._create_constants()
        GC_correction.__init__(self)
        self._deq_depth = deque(self._DEQ_LEN_ZEROS, maxlen=self._DEQ_LEN)
        self._deq_app = self._deq_depth.append
        self._deq_popl = self._deq_depth.popleft        
        self._mainlist = []
        self._mainlapp = self._mainlist.append
        self._last_ini = 0
        self.f_output = None
        self._outputlist = []
        self._fmt = '{}\t{}\t{}\t{}\t{}\t{bedcoord}\n'
        self._makeoutputfile()
        self._calltuple = namedtuple('call', 's e depth score')

    def _create_constants(self):
        self._SIZE = self.arg.SIZE
        self._OFFSET, self._NEIGHBOR = self.arg.OFFSET, self.arg.FLANKS
        self._POSITION_OFFSET = self._OFFSET+self._NEIGHBOR
        self._TOTAL_WIN_LENGTH = self._SIZE+(2*self._OFFSET)+(2*self._NEIGHBOR)
        self._CENTERINDEX = (self._SIZE-1)/2
        self._MIN_DEPTH = int(self.arg.MinDepth)
        self._DEQ_LEN = int(self.arg.DequeLength+50)
        self._DEQ_LEN_ZEROS = [0]*self._DEQ_LEN
        self._NEIGHBOR_ZEROS = [0]*self._NEIGHBOR
        self._MAX_JUMP = self._DEQ_LEN*4
        self._SPACERMEAN = float(self._NEIGHBOR+self._NEIGHBOR)

        
    def update_depth(self, record):
        if not self._last_ini:
            self._last_ini = record.pos+1
            self._last_pos = record.pos
        _jump_idx = record.pos - self._last_pos
        if _jump_idx >= self._MAX_JUMP:
            self.call_window()
            self.reset_inregion(record)
            _jump_idx = 0
        while _jump_idx:
            self._mainlapp(self._deq_popl())
            self._deq_app(0)
            _jump_idx -= 1

        self._last_pos = record.pos

        corr_depth = self._get_gc_corr_dep(record)

        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(_jump_idx, _jump_idx + count):
                    self._deq_depth[idx] += corr_depth
                _jump_idx += count
            elif cigar in (2, 3, 6):
                _jump_idx += count

    def call_window(self):
        self._mainlist.extend(self._deq_depth)
        self._call_window()

    def _check_width(self, start, end, incre):
        for idx in xrange(start, end, incre):
            val = self.window[(self._CENTERINDEX + idx)]
            if val != self.maxdepth:
                break
            self.call.append(self._CENTERINDEX+idx)

    def _depthwindows(self):
        n = len(self._mainlist)
        for idx in xrange(0, n-self._TOTAL_WIN_LENGTH+1):
            yield idx, self._mainlist[idx:(idx+self._TOTAL_WIN_LENGTH)]

    def _call_max(self, window):
        ''' docstring '''
        self.window = window
        self.maxdepth = max(self.window)
        if self.maxdepth > self._MIN_DEPTH:
            if self.window[self._CENTERINDEX] == self.maxdepth:
                self.call = [self._CENTERINDEX]
                self._check_width(1, self._CENTERINDEX, 1)
                self._check_width(-1, -self._CENTERINDEX, -1)
                # not +1 as bed files are half open at the end
                return (self.maxdepth, self._POSITION_OFFSET+min(self.call),
                        self._POSITION_OFFSET+max(self.call)+1)
        return (0, 0, 0)
        
    def _call_window(self):
        ''' docstring '''
        lasttup = ()
        for idx, win_depth in self._depthwindows():
            if win_depth[self._POSITION_OFFSET+self._CENTERINDEX] < self._MIN_DEPTH:
                continue
            if 0 in win_depth[self._POSITION_OFFSET: self._POSITION_OFFSET + self._SIZE]:
                continue
            center_depth, min_idx, max_idx = self._call_max(win_depth[self._POSITION_OFFSET:
                                                                      self._POSITION_OFFSET +
                                                                      self._SIZE])

            if center_depth:
                spacerL = sum(win_depth[:self._NEIGHBOR])
                spacerR = sum(win_depth[-self._NEIGHBOR:])

                if spacerL > self._NEIGHBOR and spacerR > self._NEIGHBOR: ## Minimum coverage of 1 in flanks
                    sizeofwindow = (max_idx-min_idx)
                    mean_spacer = (spacerL + spacerR)/self._SPACERMEAN
                    mean_spacer = mean_spacer if mean_spacer > 1 else 1
                    ## to correct for super high from the gccorrection
                    ## score = math.log(float(center_depth)/(mean_spacer*sizeofwindow))
                    score = (float(center_depth)-mean_spacer)/sizeofwindow
                    start_pos = idx+self._last_ini+min_idx
                    end_pos = idx+self._last_ini+max_idx

                    if start_pos >= self.start and end_pos <= self.end:
                        if not lasttup: # initialize
                            lasttup = self._calltuple(start_pos, end_pos, center_depth, score)
                            continue
                        # if start_pos <= lasttup.e: # the nucl dyad/center overlap
                        if start_pos <= (lasttup.e+self._SIZE-1): # the nucleosomes overlap. # use the one with highest score
                            if score > lasttup.score:
                                lasttup = self._calltuple(start_pos, end_pos, center_depth, score)
                        else:  # they do not overlap. send lasttup to outputlist
                            self._outputlist.append(lasttup)
                            lasttup = self._calltuple(start_pos, end_pos, center_depth, score)
        if lasttup:
            self._outputlist.append(lasttup)


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

    def reset_inregion(self, record):
        self._deq_depth.clear()
        self._deq_depth.extend(self._DEQ_LEN_ZEROS)
        del self._mainlist [:]
        self._mainlist.extend(self._NEIGHBOR_ZEROS)
        self._last_ini = record.pos+1-self._NEIGHBOR

    def reset_deques(self, chrom, start, end, bedcoord):
        self._deq_depth.clear()
        self._deq_depth.extend(self._DEQ_LEN_ZEROS)
        del self._mainlist [:] # everytime new bed is called. flanks are made automatically. see run func
        self._last_ini = 0
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = bedcoord
        self._outputlist = []

        
def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--MinDepth', help="..", type=int, default=5)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--GCmodel', help='..', type=str, default=None)
    parser.add_argument('--DequeLength', help="..", type=int, default=1000)
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    parser.add_argument('--MinAlignmentLength', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeSize', dest='SIZE', help="..", type=int, default=147)
    parser.add_argument('--NucleosomeFlanks', dest='FLANKS', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeOffset', dest='OFFSET', help="..", type=int, default=12)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.AlignmentFile(args.bam, "rb")
    nucl_pred_cls = Nucleosome_Prediction(args)
    flanks = int(nucl_pred_cls._TOTAL_WIN_LENGTH/2)
    for chrom, start, end, bedcoord in read_bed(args):
        nucl_pred_cls.reset_deques(chrom, start, end, bedcoord)
        start = 0 if start-flanks < 0 else start-flanks
        for record in samfile.fetch(chrom, start, end+flanks):
            if record.mapq < args.MinMappingQuality or record.is_unmapped or record.alen < args.MinAlignmentLength:
                continue  # do not analyze low quality records
            nucl_pred_cls.update_depth(record)
        nucl_pred_cls.call_window()
        nucl_pred_cls.writetofile()
    nucl_pred_cls.closefile()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
