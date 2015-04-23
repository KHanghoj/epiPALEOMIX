from __future__ import print_function
import os
import argparse
import sys
NO_BEDFILES = 3

def p_args(argv):
    parser = argparse.ArgumentParser(prog='Merge_Phasogram')
    parser.add_argument('infile', type=str)
    parser.add_argument('outdir', type=str)
    return parser.parse_known_args(argv)


def getoutnames(args):
    ''' make each bedfile filename '''
    number = 1
    fmt = '_0{}'.format
    filena, fileext = os.path.splitext(args.infile)
    while NO_BEDFILES >= number:
        yield os.path.join(args.outdir, filena + fmt(number) + fileext)
        number += 1


def getdata(args):
    ''' calculates the total nuber of lines'''
    with open(args.infile, 'r') as f_in:
        tot = sum(1 for line in f_in)
        return divmod(tot, NO_BEDFILES)


def run(args):
    ''' run produces three bedfiles for each file '''
    lineperfile, rest = getdata(args)
    with open(args.infile, 'r') as infile:
        lines = (line for line in infile)
        for f_out in getoutnames(args):
            extras = 1 if rest > 0 else 0
            with open(f_out, 'w') as outfile:
                li_per = lineperfile+extras
                # while li_per:
                #     try:
                #         outfile.write(lines.next())
                #         li_per -= 1
                #     except StopIteration:
                #         break
                while li_per:
                    outfile.write(lines.next())
                    li_per -= 1
            rest -= 1

def main(argv):
    args, _ = p_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
