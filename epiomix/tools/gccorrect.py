# !/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
#from collections import defaultdict
from epiomix.tools.commonutils import \
    Cache, \
    read_mappa
## _BUFFER = 2
_BUFFER = 0


class GCcorrect(object):
    """docstring for GCcorrect"""
    def __init__(self, arg):
        self.arg = arg
        self.samfile = pysam.AlignmentFile(self.arg.BamPath, "rb")
        self.fasta = Cache(self.arg.FastaPath)
        self.rl = self.arg.ReadLength
        self.noregions = self.arg.NoRegions
        self.reads_gc = {gc:0 for gc in xrange(self.rl+1)}
        self.reference_gc = {gc:0 for gc in xrange(self.rl+1)}
        self.halfresolution = self.arg.HalfResolution
        
    def getreads(self, chrom, start, end):
        region_size = int(end - start)
        regionseq = self.fasta.fetch_string(chrom, start-1, region_size)
        regionlen = len(regionseq)
        cov = 0
        start -= 1
        for record in self.samfile.fetch(chrom, start, end):

            if (abs(self.rl - record.alen) > self.halfresolution or  # only reads near rl
                  record.pos < start or   # not before start of region
                  record.aend-1 > end or  # not after end of region
                  record.mapq < self.arg.MinMappingQuality or
                  record.is_unmapped or
                  record.alen < self.arg.MinAlignmentLength):
                continue  # do not analyze low quality records

            cov += 1

            pos = record.aend-1-self.rl if record.is_reverse else record.pos
            pos -= start
            curr_seq = regionseq[pos:(pos+self.rl)]
            self.reads_gc[(curr_seq.count('C') +
                           curr_seq.count('G'))] += 1

        if cov:  # only allow for regions with reads
            jump = self.rl/2
            for idx in xrange(0, regionlen-self.rl+1, jump):
                curr_seq = regionseq[idx:(idx+self.rl)]
                self.reference_gc[(curr_seq.count('C') +
                                      curr_seq.count('G'))] += 1
        if cov > 100:
            self.noregions -= 1
                
    def writetofile(self):
        ''' dfs '''
        gcfmt = '{}\t{}\t{}\t{}\n'.format
        assert sum([self.reads_gc[gc] for gc in range(0, self.rl+1)]) > 2000, "Not enough reads to make the model"
        with open(self.arg.OutputFile, 'w') as f:
            for gc in range(0, self.rl+1):
                f.write(gcfmt(str(self.rl), str(gc), str(self.reads_gc[gc]),
                        str(self.reference_gc[gc])))
        self.fasta.closefile()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='GCcorrection')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('OutputFile', type=str)
    parser.add_argument('--FastaPath', type=str)
    parser.add_argument('--MappabilityPath', type=str)
    parser.add_argument('--ReadLength', help="...", type=int)
    parser.add_argument('--MappaUniqueness', help="...", type=float)
    parser.add_argument('--ChromUsed', help="...", type=str)
    parser.add_argument('--NoRegions', help="...", type=int, default=200)    
    parser.add_argument('--HalfResolution', help="...", type=int, default=4)    
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    parser.add_argument('--MinAlignmentLength', help="..", type=int, default=25)    

    return parser.parse_known_args(argv)

def run(args):
    GC = GCcorrect(args)
    mappability = args.MappaUniqueness
    last_chrom, last_end = '', -1
    runallmapparegions = False
    if args.ChromUsed.lower() == "all":
        runallmapparegions = True
    for chrom, start, end, score in read_mappa(args):
        if runallmapparegions:
            ## if "all" then entire genome analyzed or until GC.noregions hits 0
            args.ChromUsed = chrom
        if score >= mappability and args.ChromUsed == chrom and GC.noregions:
            if start-last_end < 0 and last_chrom == chrom:            # because chunks can overlap with 50%
                start += (end-start)/2
            last_chrom, last_end = chrom, end
            GC.getreads(chrom, start, end)
    GC.writetofile()

def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
