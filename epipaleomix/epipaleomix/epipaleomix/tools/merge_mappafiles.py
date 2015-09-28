from __future__ import print_function
from collections import defaultdict
import re
import argparse
import sys


def p_args(argv):
    parser = argparse.ArgumentParser(prog='Merge Mappafiles')
    parser.add_argument('analysis', type=str, help='analysisname')
    parser.add_argument('output', type=str, help='outputfile')
    parser.add_argument('infiles', nargs='+', type=str)
    return parser.parse_known_args(argv)


def unpack(gc, gc_content, readc, fastac):
    return int(gc), float(gc_content), int(readc), int(fastac)


def writetofile(out, gc, dread, dfasta):
    FMT = '{}\t{}\t{}\t{}\n'.format
    for gccont in dread.iterkeys():
        out.write(FMT(gc, gccont, dread[gccont], dfasta[gccont]))


def run(args):
    dread = defaultdict(int)
    dfasta = defaultdict(int)
    gc = 0
    with open(args.output, 'w') as f_out:
        for f in args.infiles:
            with open(f, 'r') as f_in:
                for line in f_in:
                    gc, gc_content, readc, fastac = unpack(*re.split(r'\s', line.rstrip()))
                    dread[gc_content] += readc
                    dfasta[gc_content] += fastac
        writetofile(f_out, gc, dread, dfasta)

def main(argv):
    args, unknown = p_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
