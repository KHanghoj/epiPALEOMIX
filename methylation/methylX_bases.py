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
from itertools import chain

_PLUS_STRAND_BASES = ['CG', 'TG']
_MINUS_STRAND_BASES = ['CG', 'CA']
_BASES_CHECK = 6


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
    parser.add_argument('bed', help="a bed file format with\
                        Sequences Of Interest")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_mapmethylX_bases.txt')
    return parser.parse_args(argv)


def call_ms(idxlist, last_pos, lst_dic_lastpos,
            lst_base_forward, start, dic_top, dic_lower):
    ''' docstring '''
    for idx in idxlist:
        tempdic_minus = lst_dic_lastpos[idx].pop(last_pos, {})
        top = lst_base_forward[idx].get('T', 0)+tempdic_minus.get('A', 0)
        lower = top+lst_base_forward[idx].get('C', 0)+tempdic_minus.get('G', 0)
        dic_top[idx][last_pos-start] += top
        dic_lower[idx][last_pos-start] += lower


def call_minus_ms(idxlist, last_pos, lst_dic_lastpos,
                  start, dic_top, dic_lower):
    ''' docstring '''
    lst_base_forward_temp = create_lst_dic()
    # this is needed to fix the indexing
    for idx in idxlist:
        for key in sorted(lst_dic_lastpos[idx]):
            if last_pos > key:  # this is not necessary
                call_ms([idx], key, lst_dic_lastpos, lst_base_forward_temp,
                        start, dic_top, dic_lower)


def writetofile(idxlist, dic_top, dic_lower, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    #  this is done to make sure we have all
    #  keys present in every read position
    #  even is key value is zero
    #  advantage of defaultdict
    it_keys = ((key for key in x) for x in dic_lower if x)
    keys = chain.from_iterable(it_keys)
    keys = sorted(set(keys))
    for idx in idxlist:
        for key in keys:
            f_output.write('{}\t{}\t{}\t{}\n'.format(idx,
                           key, repr(dic_top[idx][key]),
                           repr(dic_lower[idx][key])))
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


# def get_XX_index(string, forward=True, pattern='CG'):
#     start = 0
#     while True:
#         if forward:
#             idx = string.find(pattern, start)
#             if idx < 0:
#                 break
#             yield idx
#             start = idx + 1
#         else:
#             # takes the hit from the RIGHT as it the correct thing to do first
#             # for reverse strands. CG in the end is the
#             # 0th position in read.
#             idx = string.rfind(pattern, 0)
#             if idx < 0:
#                 break
#             yield idx
#             # then it cuts of the hit
#             string = string[:idx]


def get_index_forw(string, pattern='CG'):
    start = 0
    while True:

        idx = string.find(pattern, start)
        if idx < 0:
            break
        yield idx
        start = idx + 1


def get_index_rev(string, pattern='CG'):
    while True:
        # takes the hit from the RIGHT as it the correct thing to do first
        # for reverse strands. CG in the end is the
        # 0th position in read.
        idx = string.rfind(pattern, 0)
        if idx < 0:
            break
        yield idx
        # then it cuts of the hit
        string = string[:idx]


def check_lst_dic(lst):
    ''' checks which dics are not empty in list and return the idx'''
    return [idx for idx, x in enumerate(lst) if x]


def max_lst_dic(lst):
    ''' max whichs dics are not empty '''
    return max([max(x) for x in lst if x])


def create_lst_dic(nested=False, size=_BASES_CHECK-1):
    if nested:
        return [defaultdict(lambda: defaultdict(int)) for _ in range(size)]
    else:
        return [defaultdict(int) for _ in range(size)]


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    fasta = Cache(args.fastafile)
    # MS: Consider using WITH statements
    lst_dic_lastpos = create_lst_dic(nested=True)
    lst_base_forward = create_lst_dic()
    idx_standard = range(len(lst_dic_lastpos))
    dic_top = create_lst_dic()
    dic_lower = create_lst_dic()
    chrom_last = ''
    for chrom, start, end in read_bed(args):
        last_pos = -1  # reset when starting a new bedfile line
        # this is for clearing the dics
        lst_dic_lastpos = create_lst_dic(nested=True)
        lst_base_forward = create_lst_dic()
        for record in samfile.fetch(chrom, start, end):
            read_sequence = record.seq
            if record.tid != chrom_last:  # new chromosome or first record
                chrom_last = record.tid
                last_pos = -1
                pres_chrom = samfile.getrname(record.tid)

            # Call minus strand Ms with no plus strand information
            # check if some dics are not empty and max
            if check_lst_dic(lst_dic_lastpos) and \
                    max_lst_dic(lst_dic_lastpos) < last_pos:
                call_minus_ms(idx_standard, last_pos, lst_dic_lastpos,
                              start, dic_top, dic_lower)

            if record.is_reverse:  # the minus strand
                bases = read_sequence[-_BASES_CHECK:]
                pos = record.aend-_BASES_CHECK
                fast_string = fasta.fetch_string(pres_chrom,
                                                 pos, nbases=_BASES_CHECK)
                # fast_idx = list(get_XX_index(fast_string, forward=False,
                #                 pattern='CG'))
                fast_idx = list(get_index_rev(fast_string))
                if fast_idx:
                    read_idx = [x for x in fast_idx if
                                bases[x:x+2] in _MINUS_STRAND_BASES]
                    if read_idx:
                        max_pos = _BASES_CHECK - min(read_idx)
                        # because we look at the two following bases.
                        cigar_op, cigar_len = record.cigar[-1]
                        if (cigar_op == 0) and (cigar_len >= max_pos):
                            for base_idx in read_idx:
                                lst_dic_lastpos[_BASES_CHECK-base_idx-2][record.aend-2][bases[base_idx+1]] += 1
            else:  # this is for the forward strand
                bases = read_sequence[:_BASES_CHECK]
                pos = record.pos
                fast_string = fasta.fetch_string(pres_chrom,
                                                 pos, nbases=_BASES_CHECK)
                # fast_idx = list(get_XX_index(fast_string, pattern='CG'))
                fast_idx = list(get_index_forw(fast_string))
                if fast_idx:
                    read_idx = [x for x in fast_idx if
                                bases[x:x+2] in _PLUS_STRAND_BASES]
                    if read_idx:
                        max_pos = max(read_idx)+2
                        # because we look at the two following bases.
                        cigar_op, cigar_len = record.cigar[0]
                        if cigar_op == 0 and cigar_len >= max_pos:
                            if record.pos != last_pos and last_pos != -1:
                                call_ms(idx_standard, last_pos, lst_dic_lastpos,
                                        lst_base_forward, start, dic_top,
                                        dic_lower)
                                lst_base_forward = create_lst_dic()  # reset dic
                            for base_idx in read_idx:
                                lst_base_forward[base_idx][bases[base_idx]] += 1
                            last_pos = record.pos

        if check_lst_dic(lst_base_forward):  # if dics contain anything
            call_ms(idx_standard, last_pos, lst_dic_lastpos,
                    lst_base_forward, start, dic_top, dic_lower)
        if check_lst_dic(lst_dic_lastpos):
            call_minus_ms(idx_standard, last_pos, lst_dic_lastpos,
                          start, dic_top, dic_lower)
    writetofile(idx_standard, dic_top, dic_lower, args.out)
    samfile.close()
    fasta.closefile()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
