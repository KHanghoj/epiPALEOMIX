#!/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
import gzip
from collections import deque
from os.path import exists, splitext
from shutil import move
from epipaleomix.tools.commonutils import read_bed, \
    GC_correction


class Write_Depth(GC_correction):
    """docstring for Write_Depth"""
    def __init__(self, arg):
        self.arg = arg
        self._SIZE = self.arg.SIZE
        self._OFFSET, self._NEIGHBOR = self.arg.OFFSET, self.arg.FLANKS
        self._POSITION_OFFSET = self._OFFSET+self._NEIGHBOR
        self._TOTAL_WIN_LENGTH = self._SIZE+(2*self._OFFSET)+(2*self._NEIGHBOR)
        self._HALFWINDOW =  (self._TOTAL_WIN_LENGTH-1)/2  # 110
        self._CENTERINDEX = (self._SIZE-1)/2
        self._SPACERMEAN = float(self._NEIGHBOR+self._NEIGHBOR)
        GC_correction.__init__(self)
        self._f_output, self._model = None, None
        self._seq_len = int(self._TOTAL_WIN_LENGTH*2)  # maybe times 4 if small nucl
        self._zeros = [0]*self._seq_len
        self._last_pos = -self._seq_len
        self._deq_depth = deque(maxlen=self._seq_len)
        self._mainlist = list()
        self._fmt = '{}\t{}\t{}\t{}\t{}\n'
        self._makeoutputfile()

    def update_depth(self, record):
        if not self._last_ini:
            self._last_ini = record.pos+1
            self._last_pos = record.pos
        _jump_idx = record.pos - self._last_pos
        if _jump_idx:
            if _jump_idx >= self._seq_len:
                self._mainlist.extend(self._deq_depth)
                self._call_depth_scores()
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
            if self._mainlist[idx+self._HALFWINDOW]:  ## check that center is not zero depth
                yield idx, self._mainlist[idx:(idx+self._TOTAL_WIN_LENGTH)]
                
    def calcscore(self, window):
        center = window[self._HALFWINDOW]
        spacer = (sum(window[:self._NEIGHBOR])+sum(window[-self._NEIGHBOR:]))/self._SPACERMEAN
        return (center, center-spacer)
                
    def _call_depth_scores(self):
        for idx, window in self._depthwindows():
            # get to the center of the window
            pos = self._last_ini+idx+self._HALFWINDOW
            if pos >= self.start and pos <= self.end:
                depth, score = self.calcscore(window)
                self.f_output.write(self._fmt.format(self.chrom, pos,
                                                     depth, score,
                                                     self.bedcoord))

    def call_final_depth_scores(self):
        self._mainlist.extend(self._deq_depth)
        self._call_depth_scores()

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
 
    def reset_deques(self, chrom, start, end, bedcoord):
        self.chrom, self.start, self.end = chrom, start, end
        self._deq_depth = deque(self._zeros, maxlen=self._seq_len)
        self._last_ini = None
        self._mainlist = []
        self.bedcoord = bedcoord

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
    parser.add_argument('--DequeLength', help="...", type=int, default=1000)
    parser.add_argument('--NucleosomeSize', dest='SIZE', help="..", type=int, default=147)
    parser.add_argument('--NucleosomeFlanks', dest='FLANKS', help="..", type=int, default=25)
    parser.add_argument('--NucleosomeOffset', dest='OFFSET', help="..", type=int, default=12)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.Samfile(args.bam, "rb")
    Corr_Depth = Write_Depth(args)
    flanks = (Corr_Depth._TOTAL_WIN_LENGTH/2)+1
    for chrom, start, end, bedcoord in read_bed(args):
        Corr_Depth.reset_deques(chrom, start, end, bedcoord)
        start = 0 if start-flanks < 0 else start-flanks
        end += flanks
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < args.MinMappingQuality or record.is_unmapped:
                continue  # do not analyze low quality records
            Corr_Depth.update_depth(record)
        Corr_Depth.call_final_depth_scores()
    Corr_Depth.closefile()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
