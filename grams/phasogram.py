#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict

_MAX_SIZE = 3000
_MINMAPQUALI = 25
_MIN_COVERAGE = 3


class Phasogram(object):
    """docstring for Phasogram"""
    def __init__(self, arg):
        self.arg = arg
        self.outputdic = defaultdict(int)
        self.forward_dic = {}
        self.reverse_dic = {}
        self.last_tid = None
        # self.forward_dic = OrderedDict()
        # self.reverse_dic = OrderedDict()

    def _call_output(self, dic, max_lst_range=_MAX_SIZE, max_size=_MAX_SIZE):
        if dic:
            mdic = max(dic)
            while mdic - min(dic) > max_lst_range:
                old_pos = min(dic)
                old_count = dic.pop(old_pos, None)
                for current in sorted(dic.iterkeys()):
                    length = current - old_pos
                    if length >= max_size:
                        break
                        # do no know if break or contiune
                    if old_count >= _MIN_COVERAGE:
                        self.outputdic[length] += 1

    def update(self, record):
        if record.tid != self.last_tid and self.last_tid:
            self.call()
            self.reset()
        if record.is_reverse:
            curr_pos, temp_dic = record.aend, self.reverse_dic
        else:
            curr_pos, temp_dic = record.pos, self.forward_dic
        try:
            temp_dic[curr_pos] += 1
        except KeyError:
            temp_dic[curr_pos] = 1
            self._call_output(temp_dic)
        self.last_tid = record.tid

    def writetofile(self):
        with open(self.arg.out, 'w') as f_output:
            for key, value in self.outputdic.iteritems():
                f_output.write('{}\t{}\n'.format(key, value))

    def call(self):
        self._call_output(self.forward_dic, max_lst_range=0, max_size=0)
        self._call_output(self.reverse_dic, max_lst_range=0, max_size=0)

    def reset(self):
        self.forward_dic = {}
        self.reverse_dic = {}
        self.last_tid = None


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--bed', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_phasogram.txt')
    return parser.parse_args(argv)


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


def main(argv):
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    phaso = Phasogram(args)
    for chrom, start, end in read_bed(args):
        phaso.reset()
        for record in samfile.fetch(chrom, start, end):
            if record.mapq < _MINMAPQUALI:
                continue  # do not analyze low quality records
            phaso.update(record)
        phaso.call()
    phaso.writetofile()
    samfile.close()
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
