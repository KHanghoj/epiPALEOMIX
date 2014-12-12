# !/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
        # Profiling af python script:
        # $ python -m cProfile -s cumulative script.py
'''

from __future__ import print_function
import sys
import pysam
import argparse
from random import sample
from collections import defaultdict

_READ_LENGTH = 180
# _READ_LENGTH = 10
_BUFFER = 2
_N_RANDOM = 150  # 300 is maximum when frag_length is 130 and blocks are 40000
_MAPPABILITY = 0.9
# no longer random, now every position in the chunk with high mappability.


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


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('mappability_bed', help="...")
    parser.add_argument('--readlengths', help="...", default=None)
    parser.add_argument('--uniqueness', help="...", type=float, default=0.9)
    parser.add_argument('--onelength', help="...", type=int, default=31)
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_gc_content.txt')
    return parser.parse_args(argv)


def writetofile(dic_f_gc, dic_n_gc, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    for length in sorted(dic_n_gc):
        for key in sorted(dic_n_gc[length]):
            f_output.write('{}\t{}\t{}\t{}\n'.format(length, key,
                           repr(dic_f_gc[length][key]),
                           repr(dic_n_gc[length][key])))
    f_output.close()


def rand_parts(seq, n, l):
    ''' seq_fasta, no of occurences, length of occurences '''
    indices = xrange(len(seq) - (l - 1) * n)
    offset = 0
    it = iter(sorted(sample(indices, n)))
    # it = iter(indices)
    while True:
        try:
            idx = offset+it.next()
            yield idx+_BUFFER, seq[idx+_BUFFER: idx-_BUFFER+l]
            offset += l - 1
        except StopIteration:
            break


def it_fasta_seq(seq, l):
    ''' seq_fasta, length of read, return location and length'''
    indices = iter(xrange(len(seq) - (l - 1)))
    while True:
        try:
            idx = indices.next()
            yield idx+_BUFFER, seq[idx+_BUFFER: idx-_BUFFER+l]
        except StopIteration:
            break


def read_lengths(args):
    if args.readlengths:
        with open(args.readlengths, 'r') as myfile:
            for length in myfile.readlines():
                length = int(length)
                if length < 31:
                    # do not take low readlength into account
                    # should be user-defined later
                    continue
                else:
                    yield int(length)
    else:
        yield args.onelength


def read_bed(args):
    if args.mappability_bed:
        with open(args.mappability_bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')
                chrom = input_line.pop(0)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                score = float(input_line[-1])
                yield (chrom, start, end, score)
    else:
        yield (args.chrom, args.start, args.end)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = Cache(args.fastafile)
    # MS: Consider using WITH statements
    mappability = args.uniqueness
    dic_n_gc = defaultdict(lambda: defaultdict(int))
    dic_f_gc = defaultdict(lambda: defaultdict(int))
    for read_length in read_lengths(args):
        for chrom, start, end, score in read_bed(args):
            if score < mappability:
                continue
            # seq = fasta.fetch(chrom, start, end)
            seq = fasta.fetch_string(chrom, start, nbases=end-start)
            # fasta.fetch_string(chrom, start, end-start)
            for relative_pos, seq_sample in it_fasta_seq(seq, read_length):
                # this is a sliding window to create
                # the other option is random sampling.
                # see functions above, rand_parts()
                curr_start = relative_pos+start
                gc = (seq_sample.count('C')+seq_sample.count('G'))
                dic_n_gc[read_length][gc] += 1
                for record in samfile.fetch(chrom, curr_start-1,
                                            curr_start+1):
                    if record.is_reverse:
                        continue
                        ## this doesn't make sense as it's the other direction &
                        ## not within the current window in whuch GC is counted.
                        # if record.aend == curr_start:
                        #     dic_f_gc[read_length][gc] += 1
                    else:
                        if record.pos+1 == curr_start:
                            dic_f_gc[read_length][gc] += 1
    writetofile(dic_f_gc, dic_n_gc, args.out)
    samfile.close()
    fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))