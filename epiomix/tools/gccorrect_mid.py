from __future__ import print_function
import re
import sys
import argparse
from itertools import izip


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='GCcorrectionMid')
    parser.add_argument('out', type=str)
    parser.add_argument('inputs', help="..", nargs='+')
    return parser.parse_known_args(argv)


def get_data(f):
    lst_read, lst_fasta = [], []
    with open(f, 'r') as fout:
        for line in fout:
            v, p, read, fasta = re.split(r'\s', line.rstrip())
            lst_read.append(float(read))
            lst_fasta.append(float(fasta)+1)
        return v, lst_read, lst_fasta


def calc(f):
    readlength, read, fasta = get_data(f)
    sumfasta = sum(fasta)
    grandmean = sum(read)/sumfasta
    scores = []
    if not grandmean:
        return 0, 0
    for f, r in izip(fasta, read):
        scores.append(abs((r / f) - grandmean) * (f / sumfasta))
    return 0.5*sum(scores)/grandmean, readlength


def main(argv):
    args, unknown = parse_args(argv)
    valmax, optlength = 0, 0
    fmt = '{}\n'.format
    for fname in args.inputs:
        val, length = calc(fname)
        if val > valmax:
            valmax, optlength = val, length
    assert optlength > 0, "It seems there is no data available to find the optimal GCcorrection window length"
    with open(args.out, 'w') as f:
        f.write(fmt(optlength))

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
