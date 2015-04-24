from __future__ import print_function
from collections import defaultdict
import re
import argparse
import sys
# python mergephasogram.py outputnucleo test.txt test1.txt test2.txt


def p_args(argv):
    parser = argparse.ArgumentParser(prog='Merge_Phasogram')
    parser.add_argument('output', type=str, help='outputfile')
    parser.add_argument('infiles', nargs='+', type=str)
    parser.add_argument('--merge', action='store_true',
                        help='if merged output, default concatenate')
    return parser.parse_known_args(argv)


def merge_func(d, length, count):
    d[int(length)] += int(count)


def run_iter(input_gen):
    while True:
        try:
            input_gen.next()
        except StopIteration:
            break
    
    
def run(args):
    with open(args.output, 'w') as f_out:
        dic = defaultdict(int)
        for f in args.infiles:
            with open(f, 'r') as f_in:
                if args.merge:
                    # for line in f_in:
                    #     merge_func(dic, *re.split(r'\s', line.rstrip()))
                    gen = (merge_func(dic, *re.split(r'\s', line.rstrip()))
                           for line in f_in)
                else:
                    # for line in f_in:
                    #    f_out.write(line)
                    gen = (f_out.write(line) for line in f_in)
                run_iter(gen)
        if args.merge:
            for key, value in sorted(dic.iteritems()):
                f_out.write('{}\t{}\n'.format(key, value))


def main(argv):
    args, unknown = p_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
