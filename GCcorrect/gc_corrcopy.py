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

_TOTAL_FRAG_LENGTH = 140
# _TOTAL_FRAG_LENGTH = 10
_BUFFER = 2
_N_RANDOM = 100  # 300 is maximum when frag_length is 130 and blocks are 40000


# class fasta_cache(object):
#     ''' class doc '''

#     def __init__(self, filename, seq_len=40000):
#         self._samfile = pysam.Samfile(filename)
#         self._seq_len = int(seq_len)
#         self._last_chrom = None
#         self._pileups = None
#         self._last_start = None
#         self._actual_pos = None
#         self._end = None

#     def fetch_string(self, chrom, start):
#         ''' docstring '''
#         if self._last_chrom != chrom or abs(start-self._last_start) >= \
#                 self._seq_len or start >= self._end-_BUFFER:

#             self._end = start + self._seq_len
#             self._pileups = self._samfile.pileup(chrom,
#                                                 start, self._end, truncate=True)
#             self._last_start = start
#             self._last_chrom = chrom
#         self._actual_pos = start-self._last_start

#         return self._fasta_str[self._actual_pos:self._actual_pos+nbases]

def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('fastafile', help="fastafile")
    parser.add_argument('bam', help="...")
    parser.add_argument('bed', help="...")
    parser.add_argument('--unique', help="...", type=float, default=0.9)
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_gc_content.txt')
    return parser.parse_args(argv)


def writetofile(dic_f_gc, dic_n_gc, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    for key in sorted(dic_n_gc.keys()):
        f_output.write('{}\t{}\t{}\n'.format(key, repr(dic_f_gc[key]),
                       repr(dic_n_gc[key])))
    f_output.close()


# def rand_parts(seq, n, l):
#     indices = xrange(len(seq) - (l - 1) * n)
#     result = []
#     offset = 0
#     for i in sorted(random.sample(indices, n)):
#         i += offset
#         result.append(seq[i:i+l])
#         offset += l - 1
#     return result


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

    # it = read_bed(args)
    # while True:  # while true needs to have a break within
                   # or if in a function, it needs a break
                   # outside the function if in a loop
                   # like below
    #     try:
    #         chrom, start, end = it.next()

    #         ....
    #     except StopIteration:
    #         break

    # def infinite():
    #     x = 0
    #     while True:
    #         yield x
    #         x += 1


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
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
    fasta = pysam.Fastafile(args.fastafile)
    # MS: Consider using WITH statements
    dic_n_gc = defaultdict(int)
    dic_f_gc = defaultdict(int)

    for chrom, start, end, score in read_bed(args):
        if score < 0.7:
            # print(score)
            continue
        seq = fasta.fetch(chrom, start, end)
        for genomicpos, seq_sample in rand_parts(seq, _N_RANDOM,
                                                 _TOTAL_FRAG_LENGTH):
            # print(seq_sample)
            # print(genomicpos+start)
            curr_start = genomicpos+start
            gc = (seq_sample.count('C')+seq_sample.count('G'))
            dic_n_gc[gc] += 1
            # for pileupcolumn in samfile.pileup(chrom, curr_start,
            #                                    curr_start+1, truncate=True):

            #     for pileupread in pileupcolumn.pileups:
            #         if pileupread.alignment.is_reverse:  # check strand read
            #             # print(abs(pileupread.alignment.aend-1 - curr_start))
            #             if (pileupread.alignment.aend-1 == curr_start):
            #                 dic_f_gc[gc] += 1
            #                 # print('hej_reverse', curr_start, gc, seq_sample)
            #         else:
            #             forw_5_prime = (pileupcolumn.pos +
            #                             pileupread.query_position)

            #             if forw_5_prime == curr_start:
            #                 dic_f_gc[gc] += 1
            for record in samfile.fetch(chrom, curr_start-1, curr_start+1):
                if record.is_reverse:
                    if record.aend == curr_start:
                        dic_f_gc[gc] += 1
                        # print('hej_reverse', curr_start, gc, seq_sample)
                else:
                    if record.pos+1 == curr_start:
                        dic_f_gc[gc] += 1
                        # print('hej', curr_start, gc, seq_sample)
    writetofile(dic_f_gc, dic_n_gc, args.out)
    samfile.close()
    fasta.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
