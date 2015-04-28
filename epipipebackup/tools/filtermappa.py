import argparse
import sys


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help='..', type=str)
    parser.add_argument('mappa', help="..", type=float)
    return parser.parse_known_args(argv)


def main(argv):
    args, unknown = parse_args(argv)
    with open(args.inputfile, 'r') as infile:
        fmt = '{}{}'.format
        for line in infile:
            inputline = line.rstrip('\n').split('\t')
            if float(inputline[-1]) >= args.mappa:
                if 'chr' not in inputline[0]:
                    inputline[0] = fmt('chr', inputline[0])
                sys.stdout.write("\t".join(inputline)+'\n')

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
