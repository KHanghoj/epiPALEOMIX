#!/opt/local/bin/python
from __future__ import print_function
import sys
import pysam
import argparse
import gzip
from collections import deque
from os.path import exists, splitext
from shutil import move
from epiomix_commonutils import read_bed_W, \
    read_bed_WO, strtobool, Cache, GC_correction


class Write_Depth(GC_correction):
    """docstring for Write_Depth"""
    def __init__(self, arg):
        self.arg = arg
        self._fasta = Cache(self.arg.FastaPath)
        GC_correction.__init__(self)
        self._f_output, self._model = None, None
        self._GC_model_len = 0
        self._seq_len = int(self.arg.DequeLength)
        self._last_pos = -self._seq_len
        self._corrected_depth = deque(maxlen=self._seq_len)
        self._genomic_positions = deque(maxlen=self._seq_len)
        self._cache_list = list()
        self._cache_list_app = self._cache_list.append
        self._output_size = self._seq_len * 10
        self.counter = 0
        self._output_fmt = '{0}\t{1}\t{2}\t{3}\n'
        self._makeoutputfile()
        self._GCmodel_ini()

    def _check_fasta_chr(self, chrom):
        if not self.arg.FastaChromType:
            chrom = chrom.replace('chr', '')
        return chrom

    def update_depth(self, record):
        self.jump = (record.pos - self._last_pos)
        if self.jump > 0:
            if self.jump > self._seq_len:
                self.jump = self._seq_len
                pos = record.pos
            else:
                pos = self._genomic_positions[-1]+1
            if self._genomic_positions:  # retrieving depths
                self._retrieve_depth()
            self._corrected_depth.extend([0]*self.jump)
            self._genomic_positions.extend(xrange(pos, pos+self.jump))
            if self.counter > self._output_size:
                self._write_to_file()
        self._last_pos = record.pos

        if record.is_reverse:
            corr_depth = self._get_gc_corr_dep(record.aend-self._GC_model_len)
        else:
            corr_depth = self._get_gc_corr_dep(record.pos)

        deque_idx = 0
        for (cigar, count) in record.cigar:
            if cigar in (0, 7, 8):
                for idx in xrange(deque_idx, deque_idx + count):
                    self._corrected_depth[idx] += corr_depth
                deque_idx += count
            elif cigar in (2, 3, 6):
                deque_idx += count

    def _retrieve_depth(self):
        for _ in xrange(self.jump):
            dep = self._corrected_depth.popleft()
            pos = self._genomic_positions.popleft()+1
            if dep and pos >= self.start and pos <= self.end:
                self._cache_list_app((pos, dep))
                self.counter += 1

    def _write_to_file(self):
        ''' dfs '''
        # p=pos, d=depth
        for p, d in iter(self._cache_list):
            self.f_output.write(
                self._output_fmt.format(self.chrom, p,
                                        d, self.bedcoord))
        del self._cache_list[:]
        self.counter = 0

    def call_final_depths(self):
        if self._genomic_positions:
            self._retrieve_depth()
        if self._cache_list:
            self._write_to_file()

    def _makeoutputfile(self):
        ''' want to write to file every chrom, to keep scalablility'''
        if exists(self.arg.outputfile):
            pathname, extens = splitext(self.arg.outputfile)
            move(self.arg.outputfile, pathname+'_old'+extens)
        self.f_output = gzip.open(self.arg.outputfile, 'ab')
        header = '#chrom\tgenomicpos\tdepth\tbedcoord\n'
        self.f_output.write(header)

    def reset_deques(self, chrom, start, end):
        self.start, self.end = start, end
        self.chrom = self._check_fasta_chr(chrom)
        self._corrected_depth.clear()
        self._genomic_positions.clear()
        self._last_pos = -self._seq_len
        self.bedcoord = '{}_{}_{}'.format(self.chrom, self.start, self.end)

    def closefile(self):
        self.f_output.close()


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...", type=str)
    parser.add_argument('bed', help="...", type=str)
    parser.add_argument('outputfile', help='...', type=str)
    parser.add_argument('--FastaPath', help="FastaPath", type=str)
    parser.add_argument('--GCmodel', help='...', type=str, default=None)
    parser.add_argument('--FastaChromType', dest='FastaChromType')
    parser.add_argument('--BamChromType', dest='BamChromType')
    parser.add_argument('--MinMappingQuality', help="...", type=int,
                        default=25)
    parser.add_argument('--DequeLength', help="...", type=int)
    return parser.parse_known_args(argv)


def run(args):
    read_bed = read_bed_W if strtobool(args.BamChromType) else read_bed_WO
    args.FastaChromType = strtobool(args.FastaChromType)
    samfile = pysam.Samfile(args.bam, "rb")
    corrected_depth = Write_Depth(args)
    for chrom, start, end in read_bed(args):
        corrected_depth.reset_deques(chrom, start, end)
        for record in samfile.fetch(chrom, start, end):
            corrected_depth.update_depth(record)
        corrected_depth.call_final_depths()
    corrected_depth.closefile()
    return 0


def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
