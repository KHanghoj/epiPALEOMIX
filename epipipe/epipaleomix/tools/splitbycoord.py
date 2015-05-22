from __future__ import print_function
import sys
import argparse
import re

_SIZE = int(1e5)


def p_args(argv):
    parser = argparse.ArgumentParser(prog='Splitbycoord')
    parser.add_argument('inputf', type=str)
    parser.add_argument('outputf', type=str)
    return parser.parse_known_args(argv)


def unpack(c, s, e, *rest):
    return c, int(s), int(e)


def splitbycoord(args):
    ''' Splits chromosome coordinates in for every _SIZE bp
    to speed up bed split then multiprocessing'''
    bedcoord = '{}_{}_{}'.format
    fmt = '{}\t{}\t{}\t{}\n'.format
    lst = []
    lstapp = lst.append
    with open(args.inputf, 'r') as infile:
        for line in infile:
            chrom, start, end = unpack(*re.split(r'\s+', line.rstrip()))
            currbedcoord = bedcoord(chrom, start, end)
            while end - start >= _SIZE:
                lstapp((chrom, start, start + _SIZE - 1,
                        currbedcoord))
                start += _SIZE
            lstapp((chrom, start, end,
                    currbedcoord))
    return lst



def run(args):
    ''' Splits chromosome coordinates in for every _SIZE bp
    to speed up bed split then multiprocessing'''
    bedcoord = '{}_{}_{}'.format
    fmt = '{}\t{}\t{}\t{}\n'.format
    with open(args.inputf, 'r') as infile:
        with open(args.outputf, 'w') as outfile:
            for line in infile:
                chrom, start, end = unpack(*re.split(r'\s+', line.rstrip()))
                currbedcoord = bedcoord(chrom, start, end)
                while end - start >= _SIZE:
                    outfile.write(fmt(chrom, start, start + _SIZE - 1,
                                      currbedcoord))
                    start += _SIZE
                outfile.write(fmt(chrom, start, end,
                                  currbedcoord))

            
def main(argv):
    args, _ = p_args(argv)
    # run(args)
    splitbycoord(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
