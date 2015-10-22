from __future__ import print_function
from collections import defaultdict
import re
import argparse
import sys
import gzip

_HEADERS = {
    'MethylMap': '#chrom\tgenomicpos\tdeaminated\ttotal\tbedcoord\n',
    'NucleoMap': '#chrom\tstart\tend\tdepth\tscore\tbedcoord\n',
    'WriteDepth': '#chrom\tgenomicpos\tdepth\tscore\tbedcoord\n',
    'Phasogram': '#Length\tCount\n'
}
#_BUFFERSIZE = 20  # Adjust this according to how "memory efficient" you need the program to be.
_BUFFERSIZE = int(1e6)  # Adjust this according to how "memory efficient" you need the program to be.


def p_args(argv):
    parser = argparse.ArgumentParser(prog='Merge_Phasogram')
    parser.add_argument('analysis', type=str, help='analysisname')
    parser.add_argument('output', type=str, help='outputfile')
    parser.add_argument('infiles', nargs='+', type=str)
    parser.add_argument('--merge', action='store_true',
                        help='if merged output, default concatenate')
    return parser.parse_known_args(argv)

def merge_func(args):
    dic = defaultdict(int)
    for fileName in args.infiles:
        with gzip.open(fileName, 'rb') as f_in:
            for line in f_in:
                length, count = re.split(r'\s+', line.rstrip())
                dic[int(length)] += int(count)
    with gzip.open(args.output, 'wb') as f_out:
        f_out.write(_HEADERS[args.analysis])
        for key, value in sorted(dic.iteritems()):
            f_out.write('{}\t{}\n'.format(key, value))
        
def concat_func(args):
    with gzip.open(args.output, 'wb') as destFile:   ## need this prior step to add the header. cool
        destFile.write(_HEADERS[args.analysis])

    with open(args.output, 'ab') as destFile:
        for fileName in args.infiles:
            with open(fileName, 'rb') as sourceFile:
                chunk = True
                while chunk:
                    chunk = sourceFile.read(_BUFFERSIZE)
                    destFile.write(chunk)

def main(argv):
    args, unknown = p_args(argv)
    if args.merge:
        merge_func(args)
    else:
        concat_func(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
