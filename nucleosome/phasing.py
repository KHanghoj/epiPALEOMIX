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
_SIZE = 147  # the nucleosome size
_PHASING_RANGE = 300  # bases to find a new nucleosome


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bed_format', help="input is outout of nucleosome.py")
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


def writetofile(output_dic, f_name):
    ''' dfs '''
    f_output = open(f_name, 'w')
    key_values = iter(sorted((map(fun, key.split('_'))+[value] for key,
                      value in output_dic.iteritems())))
    fmt = '{0}\t{1}\t{2}\t{3}\n'
    for dat in key_values:
        f_output.write(fmt.format(*dat))
    f_output.close()


def call_score(chrom, last_start, temp_start, countscore, output_dic):
    try:
        old_start = temp_start[0]
    except IndexError:
        old_start = last_start
    key = '{}_{}_{}'.format(chrom, old_start, last_start)
    output_dic[key] = countscore


def main(argv):
    ''' docstring '''
    args = parse_args(argv)
    last_chrom = -1
    last_start = -1
    countscore = 1
    output_dic = {}
    temp_start = []
    counter_to_close = 0
    overall = 0
    last_score = -100
    for chrom, start, end, depth, score in read_bed(args):
        if abs(end - start) > 0:
            # discard all wide centers
            continue

        if chrom != last_chrom:
            last_start = 0
            countscore = 1
            temp_start = []
            last_score = -100

        shift = start - (last_start + _SIZE)
        # shift is an expression of jeg size between
        # last_start + nucleosome size minus new start
        if shift > 0:
            #  if zero or larger, nucleosome is furtheraway than
            #  147 from the former. if not. assign the new start.
            if shift <= _PHASING_RANGE:
                countscore += 1
                temp_start.append(start)
                last_start = start
                last_score = score
            else:
                call_score(chrom, last_start, temp_start,
                           countscore, output_dic)
                countscore = 1
                temp_start = []
                last_start = start
                last_score = score
        else:
            counter_to_close += 1
            if score > last_score:
                # if new score is greater than the previous assign
                # new last_start else keep the former one.
                last_start = start
                last_score = score
        overall += 1
        # last_start = start
        last_chrom = chrom
    call_score(chrom, last_start, temp_start,
               countscore, output_dic)
    writetofile(output_dic, args.out)
    print('{} called nucleosomes analyzed\n{} nucleosomes to close positioned'
          .format(overall,
          counter_to_close))
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
