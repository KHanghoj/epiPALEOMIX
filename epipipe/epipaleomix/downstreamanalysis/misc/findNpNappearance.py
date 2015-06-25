# !/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
'''
from __future__ import print_function
import sys
import pysam
import argparse
import re
import gzip
from os.path import exists
from collections import defaultdict, namedtuple
_PLUS_STRAND_BASES = ['CG', 'TG', 'GG', 'AG']
_MINUS_STRAND_BASES = ['CG', 'CA', 'CT' ,'CC']
### OBS : ACTUAL
CONV = {'C': 'G',
        'T': 'A',
        'G': 'C',
        'A': 'T'}

def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as bedfile:
            for line in bedfile:
                chrom, start, end = unpack(*re.split(r'\s+', line.rstrip()))
                yield str(chrom), int(start), int(end)

def unpack(c,s,e,*rest):
    return str(c), int(s), int(e)

                
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
            # fetch is 0-based. if 20000000 to 20000001. you get the 20000001th base
            self._fasta_str = self._fasta.fetch(chrom,
                                                start=start, end=self._end)
            self._last_start = start
            self._last_chrom = chrom
        self._actual_pos = start-self._last_start
        return self._fasta_str[self._actual_pos:self._actual_pos+nbases]

    def closefile(self):
        ''' docstring '''
        return self._fasta.close()



class FindNpN(object):
    """docstring for Methyl_Level"""
    def __init__(self, arg):
        self.arg = arg
        self._ReadBases = self.arg.ReadBases
        self._fastaCache = Cache(self.arg.FastaPath)
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.pat = re.compile('CG')
        self.rows = list()
        self.rowsapp = self.rows.append
        self.na_tup = namedtuple('row', ('pos A C G T'))
        self._makeoutputfile()
        self._row_size = int(1e6)
        self.forw_five = {'inbases': _PLUS_STRAND_BASES,
                          'basepos': 0,
                          'skip': self.arg.SkipFivePrime}
        self.rev_five = {'inbases': _MINUS_STRAND_BASES,
                         'basepos': 1,
                         'skip': self.arg.SkipFivePrime}

    def _updatefunc(self):
        if self.record.is_reverse:
            self._rightpart(**self.rev_five)
        else:
            self._leftpart(**self.forw_five)


    def reset_dict(self, chrom, start, end):
        self.dic_pos = defaultdict(lambda: defaultdict(int))
        self.start, self.end, self.chrom = start, end, chrom
        del self.rows[:]  # clear all data
        self.last_end = 0

    def _getindexes(self, bases_str):
        ''' returns the 0-based indeces of fasta read'''
        return (m.start() for m in self.pat.finditer(bases_str))

    def update(self, record):
        self.record = record
        if self.last_end < self.record.pos:
            self._call_ms()
            self.dic_pos.clear()
            if len(self.rows) > self._row_size:
                self._writetofile()
        self._updatefunc()
        self.last_end = self.record.aend

    def _call_ms(self):
        for pos, basescore in sorted(self.dic_pos.iteritems()):
            if pos >= self.start and pos <= self.end:  # make sure overlapping genomic sites do not get counted twice
                self.rowsapp(self.na_tup(pos,
                                         basescore.get('A', 0),
                                         basescore.get('C', 0),
                                         basescore.get('G', 0),
                                         basescore.get('T', 0)))

    def call_final_ms(self):
        self._call_ms()
        self._writetofile()

    def _prep(self, curr_pos, skip=0):
        ## skip is only used in right side functions
        # skip is 1-based always
        fast_string = self._fastaCache.fetch_string(self.chrom, curr_pos,
                                               self._ReadBases-skip)
        for fast_idx in self._getindexes(fast_string):
            yield fast_idx, self._ReadBases - fast_idx

    def _rightpart(self, inbases, basepos, skip):
        curr_pos = self.record.aend-self._ReadBases
        cigar_op, cigar_len = self.record.cigar[-1]
        bases = self.record.seq[-self._ReadBases:]
        for fast_idx, inverse_idx in self._prep(curr_pos, skip):
            if (cigar_op == 0 and cigar_len >= inverse_idx and
                    bases[fast_idx:fast_idx+2] in inbases):
                self.dic_pos[curr_pos+fast_idx+1][CONV[bases[fast_idx+basepos]]] += 1
        
    def _leftpart(self, inbases, basepos, skip):
        curr_pos = self.record.pos
        cigar_op, cigar_len = self.record.cigar[0]
        bases = self.record.seq[:self._ReadBases]
        for fast_idx, _ in self._prep(curr_pos):
            if (fast_idx >= skip and cigar_op == 0 and
                    cigar_len >= fast_idx and
                    bases[fast_idx:fast_idx+2] in inbases):
                self.dic_pos[curr_pos+fast_idx+1][bases[fast_idx+basepos]] += 1            

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalable'''
        assert not exists(self.arg.outputfile), 'File exists already. Delete or move "%s" to run analyses with same outputname' % self.arg.outputfile
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
        self.f_output.write('#chr\tpos\tA\tC\tG\tT\n')
        self.fmt = "{chr}\t{r.pos}\t{r.A}\t{r.C}\t{r.G}\t{r.T}\n".format

    def _writetofile(self):
        for row in self.rows:
            self.f_output.write(self.fmt(r=row, chr=self.chrom))
        del self.rows[:]  # EMPTY LIST AFTER WRITING TO FILE

    def closefiles(self):
        self._fastaCache.closefile()
        self.f_output.close()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('FastaPath', help="FastaPath", type=str)
    parser.add_argument('bam', help="..", type=str)
    parser.add_argument('bed', help="..", type=str)
    parser.add_argument('outputfile', help='..', type=str)
    parser.add_argument('--ReadBases', help="..", type=int, default=2)
    parser.add_argument('--MinMappingQuality', help="..", type=int, default=25)
    parser.add_argument('--SkipFivePrime', help="..", type=int, default=0)
    return parser.parse_known_args(argv)


def run(args):
    samfile = pysam.Samfile(args.bam, "rb")
    Met_Score =  FindNpN(args)
    for chrom, start, end in read_bed(args):
        Met_Score.reset_dict(chrom, start, end)
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < args.MinMappingQuality:
                continue  # do not analyze low quality records
            Met_Score.update(record)
        Met_Score.call_final_ms()
    Met_Score.closefiles()
    samfile.close()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    print(args)
    print(unknown)
    run(args)
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
