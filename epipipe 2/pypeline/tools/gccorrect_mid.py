from __future__ import print_function
import re
import sys
import numpy as np
import argparse
# python tools/gccorrect_mid.py tralaa.txt
#     temp/Saqqaq/Saqqaq_GCcorrect_*


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
            lst_fasta.append(float(fasta))
        return v, np.asarray(lst_read), np.asarray(lst_fasta)+1


def calc(f):
    readlength, read, fasta = get_data(f)
    grandmean = sum(read)/sum(fasta)
    dat = sum(abs(read / fasta - grandmean) * (fasta / sum(fasta)))
    return 0.5*dat/grandmean, readlength


def main(argv):
    args, unknown = parse_args(argv)
    valmax, optlength = 0, 0
    # fmt = '{}\t{}\n'.format
    fmt = '{}\n'.format
    for fname in args.inputs:
        val, length = calc(fname)
        if val > valmax:
            valmax, optlength = val, length
    with open(args.out, 'w') as f:
        # f.write(fmt(valmax, optlength))
        f.write(fmt(optlength))

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
