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
from collections import Counter, defaultdict

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


def writetofile(dic_f_gc, dic_n_gc, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    for key in sorted(dic_n_gc.keys()):
        f_output.write('{}\t{}\t{}\n'.format(key, repr(dic_f_gc[key]),
                       repr(dic_n_gc[key])))
    f_output.close()


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    # f_output = open(args.out, 'w')
    fasta = pysam.Fastafile(args.fastafile)
    # MS: Consider using WITH statements
    chrom = None
    dic_fasta_gc = Counter()
    dic_n_gc = defaultdict(int)
    dic_f_gc = defaultdict(int)
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
                    actual_read = fasta_wind[_BEGIN_A:_END_M]  # removbegin&end
                    dic_fasta_gc.update(actual_read)
                    gc = dic_fasta_gc['G']+dic_fasta_gc['C']
                    temp_start = start + current_pos + _BEGIN_A
                    # f_output.write('{}\t{}\t{}\n'.format(chrom,
                    #                repr(temp_start), repr(gc)))

                    # fasta_wind.pop(0)
                    dic_n_gc[gc] += 1
                    # temp_end = temp_start+_TOTAL_FRAG_LENGTH
                    for record in samfile.fetch(chrom, temp_start,
                                                temp_start+1):
                        if record.pos+1 == temp_start:
                            print(record.pos+1, temp_start)
                            dic_f_gc[gc] += 1
                    print()
                    dic_fasta_gc.clear()
                    del fasta_wind[:]

    writetofile(dic_f_gc, dic_n_gc, args.out)

    # f_output.close()
    samfile.close()
    fasta.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
