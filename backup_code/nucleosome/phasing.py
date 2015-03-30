#!/opt/local/bin/python
'''  Object: call the nucleosome center if exists and return key=POSITION,
value=DEPTH. It extends the values around the centre of the nucleosome if
identical maximal depth
Example:
test.bam --chrom 22 --start 18500000 --end 18513769 --out dingdong
Thi is with and empty '\n' bed file.
nucleomap_bedformat.py test.bam test.bed --chrom 22 --start 20000000 --end 20200000 --out tra12
'''
from __future__ import print_function
import sys
import argparse

#### Constants:
_SIZE = 147  # the nucleosome size 147 base pairs
_PHASING_RANGE = 400  # bases to find a new nucleosome


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_format', help="input is outout of nucleosome.py")
    parser.add_argument('--UCSCformat', help='...', default=False)
    parser.add_argument('--out', help='...', default='out_phasing.txt')

    return parser.parse_args(argv)


def read_bed(args):
    with open(args.bed_format, 'r') as myfile:
        for line in myfile.readlines():
            input_line = (line.rstrip('\n')).split('\t')
            chrom = input_line.pop(0).replace('chr', '')
            start = int(input_line.pop(0))
            end = int(input_line.pop(0))
            depth = int(input_line.pop(0))
            score = float(input_line.pop(0))
            yield (chrom, start, end, depth, score)


def fun(item):
    try:
        return int(item)
    except ValueError:
        return str(item)


def writetofile(output_dic, f_name, UCSCformat):
    ''' dfs '''
    f_output = open(f_name, 'w')
    key_values = iter(sorted((map(fun, key.split('_'))+[value] for key,
                      value in output_dic.iteritems())))
    if UCSCformat:
        fmt = 'chr{0}:{1}-{2}\t{3}\n'
    else:
        fmt = '{0}\t{1}\t{2}\t{3}\n'
    for dat in key_values:
        f_output.write(fmt.format(*dat))
    f_output.close()


def call_score(chrom, lst_start, output_dic):
    key = '{}_{}_{}'.format(chrom, lst_start[0], lst_start[-1])
    output_dic[key] = len(lst_start)


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    last_chrom = ''
    output_dic = {}
    lst_start = []
    counter_to_close = 0
    overall = 0
    for chrom, start, end, depth, score in read_bed(args):
        if abs(end - start) > 0:
            # discard all wide centers
            continue

        if chrom != last_chrom:
            if lst_start:
                call_score(chrom, lst_start, output_dic)
            lst_start = [start]

        # shift = start - last_start
        shift = start - lst_start[-1]
        if shift >= _SIZE:
            if shift <= _PHASING_RANGE:
                lst_start.append(start)
            else:
                call_score(chrom, lst_start, output_dic)
                lst_start = [start]
        else:
            counter_to_close += 1
        overall += 1
        last_chrom = chrom
    call_score(chrom, lst_start, output_dic)
    writetofile(output_dic, args.out, args.UCSCformat)
    print('{} 0bp-wide nucleosomes analyzed\n{} nucleosomes to close positioned'
          .format(overall, counter_to_close))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
