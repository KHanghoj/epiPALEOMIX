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
from collections import defaultdict

_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']


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
    parser.add_argument('--bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_mapmethyl.txt')
    return parser.parse_args(argv)


def call_ms(last_pos, dic_lastpos, dic_base_forward, start, dic_top, dic_lower):
    ''' docstring '''
    tempdic_minus = dic_lastpos.pop(last_pos, {})
    top = dic_base_forward.get('T', 0)+tempdic_minus.get('A', 0)
    lower = top+dic_base_forward.get('C', 0)+tempdic_minus.get('G', 0)
    dic_top[last_pos-start] += top
    dic_lower[last_pos-start] += lower


def call_minus_ms(last_pos, dic_lastpos, start, dic_top, dic_lower):
    ''' docstring '''
    for key in sorted(dic_lastpos):
        if last_pos > key:  # this is not necessary
            call_ms(key, dic_lastpos, {}, start, dic_top, dic_lower)


def writetofile(dic_top, dic_lower, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    for key in sorted(dic_top):
        f_output.write('{}\t{}\t{}\n'.format(key, repr(dic_top[key]),
                       repr(dic_lower[key])))
    f_output.close()


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as myfile:
            for line in myfile.readlines():
                input_line = (line.rstrip('\n')).split('\t')[:3]
                chrom = input_line.pop(0).replace('chr', '')
                start = int(input_line.pop(0))
                end = int(input_line.pop(0))
                yield (chrom, start, end)
    else:
        yield (args.chrom, args.start, args.end)


def get_XX_index(string, pattern='CG'):
    start = 0
    while True:
        idx = string.find(pattern, start)
        if idx < 0:
            break
        yield idx
        start = idx + 1


# def check_read(string, indices, bases):
#     for idx in indices:
#         if string in bases:
#             yield True
#         else: yield False
#     it = iter(indices)
#     while True:
#         try:
#             idx = it.next
#             if not string[idx:idx+2] in _PLUS_STRAND_BASES:
#                 break
#             # yield string[idx:idx+2]
#         except StopIteration:
#             break


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = Cache(args.fastafile)
    # MS: Consider using WITH statements
    chrom = None
    # dic_lastpos = defaultdict(lambda: defaultdict(int))
    lst_dic_lastpos = [defaultdict(lambda: defaultdict(int)) for _ in range(6)]
    # dic_base_forward = defaultdict(int)
    lst_base_forward = [defaultdict(int) for _ in range(6)]
    dic_top = defaultdict(int)
    dic_lower = defaultdict(int)
    for chrom, start, end in read_bed(args):
        last_pos = -1  # reset when starting a new bedfile line
        # dic_lastpos.clear()
        # dic_base_forward.clear()
        lst_dic_lastpos = [defaultdict(lambda: defaultdict(int)) for _ in range(6)]
        lst_base_forward = [defaultdict(int) for _ in range(6)]
        for record in samfile.fetch(chrom, start, end):
            read_sequence = record.seq
            if record.tid != chrom:  # new chromosome or first record
                chrom = record.tid
                last_pos = -1
                pres_chrom = samfile.getrname(record.tid)

            # Call minus strand Ms with no plus strand information
            if dic_lastpos and max(dic_lastpos) < last_pos:
                call_minus_ms(last_pos, dic_lastpos,
                              start, dic_top, dic_lower)

            if record.is_reverse:  # the minus strand
                continue  # only TEMPORARY
                # if read_sequence[-2:] in _MINUS_STRAND_BASES:
                #     # last two bases ok
                #     pos = record.aend-2  # end of read minus 2 bases
                #     if 'CG' in fasta.fetch_string(pres_chrom, pos):
                #         cigar_op, cigar_len = record.cigar[-1]
                #         if (cigar_op == 0) and (cigar_len >= 2):
                #             dic_lastpos[record.aend-2][read_sequence[-1]] += 1

            else:  # this is for the forward strand
                bases = read_sequence[:6]
                fast_string = fasta.fetch_string(pres_chrom,
                                                 record.pos, nbases=6)
                fast_idx = list(get_XX_index(fast_string, pattern='CG'))
                read_idx = (x for x in fast_idx if
                            bases[x:x+2] in _PLUS_STRAND_BASES)

                if fast_idx and read_idx:
                    matchidx = max(read_idx)
                    cigar_op, cigar_len = record.cigar[0]
                    if cigar_op == 0 and cigar_len >= matchidx:
                        if record.pos != last_pos and last_pos != -1:
                            call_ms(last_pos, dic_lastpos,
                                    dic_base_forward, start, dic_top,
                                    dic_lower)
                            dic_base_forward.clear()
                        # for base_idx in read_idx:
                            # dic_base_forward[read_sequence[base_idx]] += 1
                        last_pos = record.pos
        if dic_base_forward:  # this checks if dic contains anything
            call_ms(last_pos, dic_lastpos, dic_base_forward, start, dic_top,
                    dic_lower)
        if dic_lastpos:
            call_minus_ms(last_pos, dic_lastpos, start, dic_top, dic_lower)
    writetofile(dic_top, dic_lower, args.out)
    samfile.close()
    fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
