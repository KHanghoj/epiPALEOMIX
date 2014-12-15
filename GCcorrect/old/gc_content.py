# !/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
 python ~/research/projects/epiomix/methylation/methyl_fetch.py chrome.fa test.bam --chrom 22 --start 18100000 --end 20100000 --out new.txt
 Rscript -e "a=read.table('new.txt') ;summary(a)"
 Rscript -e "a=read.table('new.txt');b=sum(a[,1])/sum(a[,2]);print(b)"
        # Profiling af python script:
        # $ python -m cProfile -s cumulative script.py
'''

from __future__ import print_function
import sys
import pysam
import argparse
from collections import Counter

# _PLUS_STRAND_BASES = ['CG', 'TG']
# _MINUS_STRAND_BASES = ['CG', 'CA']

_TOTAL_FRAG_LENGTH = 31
_BEGIN_A = 2
_END_M = -2
_BASES = ['A', 'C', 'G', 'T']


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


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    f_output = open(args.out, 'w')
    fasta = pysam.Fastafile(args.fastafile)
    # MS: Consider using WITH statements
    chrom = None
    dic_fasta_gc = Counter()
    fasta_wind = []
    with open(args.bed, 'r') as bedfile_f:
        for line in bedfile_f.readlines():
            input_line = (line.rstrip('\n')).split('\t')
            if float(input_line[-1]) < args.unique:
                # removes regions with low uniqueness
                continue
            del fasta_wind[:]
            dic_fasta_gc.clear()
            current_pos = -1
            # it is the users responsibility to input bed format
            # identical to BAM format
            try:
                chrom = input_line.pop(0)
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
            except (ValueError, IndexError):
                print('something si not right')
                sys.exit()
                # chrom = args.chrom
                # start = args.start
                # end = args.end
            for fasta_base in fasta.fetch(chrom, start, end):
                current_pos += 1
                if fasta_base not in _BASES:
                    continue
                fasta_wind.append(fasta_base)
                if len(fasta_wind) == _TOTAL_FRAG_LENGTH:
                    actual_read = fasta_wind[_BEGIN_A:_END_M]  # remov begin&end
                    dic_fasta_gc.update(actual_read)
                    gc = dic_fasta_gc['G']+dic_fasta_gc['C']
                    f_output.write('{}\t{}\t{}\n'.format(chrom,
                                   repr(start+current_pos), repr(gc)))
                    dic_fasta_gc.clear()
                    # fasta_wind.pop(0)
                    del fasta_wind[:]

    f_output.close()
    samfile.close()
    fasta.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
