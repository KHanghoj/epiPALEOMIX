#!/usr/bin/env python
from __future__ import print_function
import argparse
import sys
import re
_SIZE = int(1e6)
# BEDCOORD = '{}_{}_{}'.format
BEDCOORD = '{}_{}_{}_{strand}'.format

def p_args(argv):
    parser = argparse.ArgumentParser(prog='split bedfiles')
    parser.add_argument('infile', type=str)
    parser.add_argument('outfiles', nargs='+', type=str)
    return parser.parse_known_args(argv)


def unpack(c, s, e, *rest):
    if rest and re.search('^\S+_\d+_\d+_[-+]$', rest[0]):
        bedcoord = rest[0]  ## it is highly likely to be a bedcoord
    else:
        curr_strand = re.search(r"(\s+|^)(?P<strandtype>[+-])(\s+|$)"," ".join(rest))
        if curr_strand:
            bedcoord = BEDCOORD(c, s, e,
                                strand=curr_strand.group("strandtype"))
        else:
            bedcoord = BEDCOORD(c, s, e, strand='+')
    return str(c), int(s), int(e), bedcoord


def splitbycoord(args):
    ''' Splits chromosome coordinates in for every _SIZE bp
    to speed up bed split then multiprocessing'''
    fmt = '{}\t{}\t{}\t{}\n'.format
    lst = []
    lstapp = lst.append
    with open(args.infile, 'r') as infile:
        for line in infile:
            chrom, start, end, currbedcoord = unpack(*re.split(r'\s+', line.rstrip()))
            while end - start >= _SIZE:
                lstapp(fmt(chrom, start, start + _SIZE - 1,
                           currbedcoord))
                start += _SIZE
            if end-start > 0:
                lstapp(fmt(chrom, start, end,
                           currbedcoord))
    return iter(lst), len(lst)

    

def run(args):
    ''' split bedfile into subbedfiles. divmod controls the number of lines per file '''
    no_subbed = len(args.outfiles)  # noumber of sub files is based on no of outputfiles given as argument
    infiledata, nolines = splitbycoord(args)
    lineperfile, rest = divmod(nolines, no_subbed)
    for f_out in args.outfiles:
        extras = 1 if rest > 0 else 0
        with open(f_out, 'w') as outfile:
            li_per = lineperfile+extras
            while li_per:
                outfile.write(infiledata.next())
                li_per -= 1
            rest -= 1

def main(argv):
    args, _ = p_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
