import argparse
import sys


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfile', help='..', type=str)
    parser.add_argument('mappauniqueness', help="..", type=float)
    return parser.parse_known_args(argv)


def main(argv):
    args, unknown = parse_args(argv)
    with open(args.inputfile, 'r') as infile:
        count = 0
        for line in infile:
            inputline = line.rstrip("\n").split("\t")
            if float(inputline[-1]) >= args.mappauniqueness:
                count += 1
                sys.stdout.write("\t".join(inputline[:3])+"\n") ## write coordinates only
        assert count, "No Regions in the genome displayed complexity levels equal to or higher than %s\n Recommended uniqueness score 0.75 - 0.9" % (args.mappauniqueness)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
