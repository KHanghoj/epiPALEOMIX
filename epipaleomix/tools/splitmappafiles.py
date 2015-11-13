from __future__ import print_function
import argparse
import sys
import re


def p_args(argv):
    parser = argparse.ArgumentParser(prog='split mappa')
    parser.add_argument('infile', type=str)
    parser.add_argument('outfiles', nargs='+', type=str)
    return parser.parse_known_args(argv)


def unpack(c, s, e, *rest):
    return c, s, e, rest[-1]


def getdata(args):
    fmt = '{}\t{}\t{}\t{}\n'.format
    lst = []
    lstapp = lst.append
    with open(args.infile, 'r') as infile:
        for line in infile:
            chrom, start, end, score = unpack(*re.split(r'\s+', line.rstrip()))
            lstapp(fmt(chrom, start, end,
                       score))
    return iter(lst), len(lst)


def run(args):
    ''' run produces three bedfiles for each file '''
    no_subbed = len(args.outfiles)  # noumber of sub files is based on no of outputfiles given as argument
    infiledata, nolines = getdata(args)
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
