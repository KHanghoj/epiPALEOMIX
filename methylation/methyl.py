#!/opt/local/bin/python
'''  Object: To find the methylation value from a region. the methylation
score (Ms)
'''

from __future__ import print_function
import sys
import pysam
import math
import argparse

#### Constants:




def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('name', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_mapmethyl.txt')
    return parser.parse_args(argv)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    print(args.name)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
