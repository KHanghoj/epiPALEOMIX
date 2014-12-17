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
# import gzip

# _READ_LENGTH = 180
# _READ_LENGTH = 10
_BUFFER = 2
# _N_RANDOM = 150  # 300 is maximum when frag_length is 130 and blocks are 40000
# _MAPPABILITY = 0.9
# no longer random, now every position in the chunk with high mappability.


class Cache_fasta(object):
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
    parser.add_argument('--onelength', help="...", type=int, default=56)
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_gc_content.txt')
    return parser.parse_args(argv)


def writetofile(dic_f_gc, dic_n_gc, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    for length in sorted(dic_n_gc):
        # for key in sorted(dic_n_gc[length]):
        for key in range(0, length+1):
            f_output.write('{}\t{}\t{}\t{}\n'.format(length, key,
                           repr(dic_f_gc[length][key]),
                           repr(dic_n_gc[length][key])))
    f_output.close()


def rand_parts(seq, n, l):
    ''' seq_fasta, no of occurences, length of occurences '''
    indices = xrange(len(seq) - (l - 1) * n)
    offset = 0
    it = iter(sorted(sample(indices, n)))
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
                    yield length
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


def update(pos, dic):
    dic[pos] = dic.get(pos, 0) + 1


def call_five_prime_pos(relative_pos, seq_sample, start, read_length,
                        dic_plus, dic_minus, dic_n_gc, dic_f_gc):
    curr_start = relative_pos+start
    curr_end = curr_start + len(seq_sample)
    gc = seq_sample.count('C')+seq_sample.count('G')
    dic_n_gc[read_length][gc] += 2
    ## this is two as both ways directions in one go.
    plus_count = dic_plus.get(curr_start-1, None)
    minus_count = dic_minus.get(curr_end, None)
    # plus_count = dic_plus.pop(curr_start-1, None)
    # minus_count = dic_minus.pop(curr_end-1, None)
    if plus_count:
        dic_f_gc[read_length][gc] += plus_count
    if minus_count:
        dic_f_gc[read_length][gc] += minus_count


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = Cache_fasta(args.fastafile)
    # MS: Consider using WITH statements
    mappability = args.uniqueness
    dic_n_gc = defaultdict(lambda: defaultdict(int))
    dic_f_gc = defaultdict(lambda: defaultdict(int))
    last_chrom, last_end = '', -1

    for chrom, start, end, score in read_bed(args):
        if score < mappability:
            continue
        if start-last_end < 0 and last_chrom == chrom:
            start = start + ((end-start)/2)
        last_chrom = chrom
        last_end = end
        dic_plus = {}
        dic_minus = {}
        lengthofregion = end-start
        averagehit = lengthofregion*0.0002
        # have at least one read for every 5Kb.

        records = samfile.fetch(chrom, start, end)
        # here we get all postions in minus strand and plus strand in a dict
        [update(record.aend, dic_minus) if record.is_reverse else
            update(record.pos, dic_plus) for record in records]

        # if sum(dic_minus.values())+sum(dic_plus.values())\
        #         < (lengthofregion)/100:
            # need to have at least one read for every 100 bases on average
            # continue

        if len(dic_minus)+len(dic_plus) < averagehit:
            # need to have at least one read
            # starting at average of every 1KB
            continue
        print(chrom, start, end, len(dic_minus)+len(dic_plus))

        seq = fasta.fetch_string(chrom, start, nbases=lengthofregion)
        for read_length in read_lengths(args):
            # now we run through all the different read lengths
            # can be in the beginnig as well.
            # if in beginning we can pop the dict.
            for relative_pos, seq_sample in it_fasta_seq(seq, read_length):
                call_five_prime_pos(relative_pos, seq_sample,
                                    start, read_length, dic_plus, dic_minus,
                                    dic_n_gc, dic_f_gc)
    writetofile(dic_f_gc, dic_n_gc, args.out)
    samfile.close()
    fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
