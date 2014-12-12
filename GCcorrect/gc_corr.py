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

_TOTAL_FRAG_LENGTH = 180
# _TOTAL_FRAG_LENGTH = 10
_BUFFER = 2
_N_RANDOM = 150  # 300 is maximum when frag_length is 130 and blocks are 40000
_MAPPABILITY = 0.9


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
        if score < _MAPPABILITY:
            continue
        seq = fasta.fetch(chrom, start, end)
        for genomicpos, seq_sample in rand_parts(seq, _N_RANDOM,
                                                 _TOTAL_FRAG_LENGTH):
            curr_start = genomicpos+start
            gc = (seq_sample.count('C')+seq_sample.count('G'))
            dic_n_gc[gc] += 1
            for record in samfile.fetch(chrom, curr_start-1, curr_start+1):
                if record.is_reverse:
                    if record.aend == curr_start:
                        dic_f_gc[gc] += 1
                else:
                    if record.pos+1 == curr_start:
                        dic_f_gc[gc] += 1
    writetofile(dic_f_gc, dic_n_gc, args.out)
    samfile.close()
    fasta.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
