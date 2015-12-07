#!/usr/bin/env python
from __future__ import print_function
import sys
import pysam
import argparse
import gzip
from collections import deque, namedtuple
from os.path import exists, splitext
from shutil import move
from epiomix.tools.commonutils import read_bed, \
    GC_correction


class Write_Depth(GC_correction):
    """docstring for Write_Depth"""
    def __init__(self, arg):
        self.arg = arg
        self._create_constants()
        GC_correction.__init__(self)
        self._f_output = None
        self._deq_depth = deque(self._DEQ_LEN_ZEROS, maxlen=self._DEQ_LEN)
        self._deq_app = self._deq_depth.append
        self._deq_popl = self._deq_depth.popleft        
        self._mainlist = []
        self._mainlapp = self._mainlist.append
        self._outputlist = []
        self._calltuple = namedtuple('call', 'pos depth score')
        self._fmt = '{}\t{}\t{}\t{}\t{bedcoord}\n'
        self._makeoutputfile()

    def _create_constants(self):
        self._SIZE = self.arg.SIZE
        self._OFFSET, self._NEIGHBOR = self.arg.OFFSET, self.arg.FLANKS
        self._POSITION_OFFSET = self._OFFSET+self._NEIGHBOR
        self._TOTAL_WIN_LENGTH = self._SIZE+(2*self._OFFSET)+(2*self._NEIGHBOR)
        self._HALFWINDOW = int((self._TOTAL_WIN_LENGTH-1)/2)  ## +1
        self._CENTERINDEX = (self._SIZE-1)/2
        self._DEQ_LEN = int(self.arg.DequeLength+50)
        self._DEQ_LEN_ZEROS = [0]*self._DEQ_LEN
        self._HALFWINDOW_ZEROS = [0]*self._HALFWINDOW
        self._MAX_JUMP = self._DEQ_LEN*4
        self._SPACERMEAN = float(self._NEIGHBOR+self._NEIGHBOR)

    def update_depth(self, record):
        if not self._last_ini:
            self._last_ini = record.pos+1
            self._last_pos = record.pos
        _jump_idx = record.pos - self._last_pos
        if _jump_idx >= self._MAX_JUMP:
            self.call_depths()
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
                
    def call_depths(self):
        if sum(self._mainlist) or sum(self._deq_depth):
            self._mainlist.extend(self._deq_depth)
            self._mainlist.extend(self._HALFWINDOW_ZEROS)
            self._call_depth_scores()

    def _depthwindows(self):
        n = len(self._mainlist)
        for idx in xrange(0, n-self._TOTAL_WIN_LENGTH+1):
            if self._mainlist[idx+self._HALFWINDOW]:
                yield idx, self._mainlist[idx:(idx+self._TOTAL_WIN_LENGTH)]

    def _calcscore(self, window):
        center = window[self._HALFWINDOW]
        spacer = (sum(window[:self._NEIGHBOR])+sum(window[-self._NEIGHBOR:]))/self._SPACERMEAN
        # spacer = spacer if spacer > 1 else 1
        return (center, center-spacer)
                
    def _call_depth_scores(self):
        for idx, window in self._depthwindows():
            pos = self._last_ini+idx+self._HALFWINDOW
            if pos >= self.start and pos <= self.end:
                depth, score = self._calcscore(window)
                self._outputlist.append(self._calltuple(pos, depth, score))

    def writetofile(self):
        for dat in self._outputlist:
            self.f_output.write(self._fmt.format(self.chrom, *dat,
                                                 bedcoord=self.bedcoord))
        self._outputlist = [] 

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
 
    def reset_inregion(self, record):
        self._deq_depth.clear()
        self._deq_depth.extend(self._DEQ_LEN_ZEROS)
        del self._mainlist [:]
        self._mainlist.extend(self._HALFWINDOW_ZEROS)
        self._last_ini = record.pos+1-self._HALFWINDOW

    def reset_deques(self, chrom, start, end, bedcoord):
        self._deq_depth.clear()
        self._deq_depth.extend(self._DEQ_LEN_ZEROS)
        del self._mainlist [:] # everytime new bed is called. flanks are made automatically. see run func
        self._last_ini = 0
        self.start, self.end, self.chrom = start, end, chrom
        self.bedcoord = bedcoord
        self._outputlist = []

        
    def closefile(self):
        self.f_output.close()
        try:
            self._fasta_dat.closefile()
        except AttributeError:
            pass


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...", type=str)
    parser.add_argument('bed', help="...", type=str)
    parser.add_argument('outputfile', help='...', type=str)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--GCmodel', help='...', type=str, default=None)
    parser.add_argument('--MinMappingQuality', help="...", type=int, default=25)
    parser.add_argument('--MinAlignmentLength', help="...", type=int, default=25)
    parser.add_argument('--DequeLength', help="...", type=int, default=1000)
    parser.add_argument('--NucleosomeSize', dest='SIZE', help="..", type=int, default=147)
    parser.add_argument('--NucleosomeFlanks', dest='FLANKS', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeOffset', dest='OFFSET', help="..", type=int, default=12)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.AlignmentFile(args.bam, "rb")
    Corr_Depth = Write_Depth(args)
    # flanks = (Corr_Depth._TOTAL_WIN_LENGTH/2)+1
    flanks = Corr_Depth._TOTAL_WIN_LENGTH
    for chrom, start, end, bedcoord in read_bed(args):
        Corr_Depth.reset_deques(chrom, start, end, bedcoord)
        start = 0 if start-flanks < 0 else start-flanks
        end += flanks
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < args.MinMappingQuality or record.is_unmapped or record.alen < args.MinAlignmentLength:
                continue  # do not analyze low quality records
            Corr_Depth.update_depth(record)
        Corr_Depth.call_depths()
        Corr_Depth.writetofile()
    Corr_Depth.closefile()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
