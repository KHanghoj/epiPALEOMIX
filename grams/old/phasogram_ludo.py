#!/opt/local/bin/python
''' Object: To calculate the phasogram between 5' end
of mapped reads's start positions
aligning in same orientation within a 1000bp window.
'''

from __future__ import print_function
import sys
import pysam
import argparse
from collections import defaultdict

_MIN_DEPTH = 5
_MAX_SIZE = 1000
_MINMAPQUALI = 30


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser()
    parser.add_argument('bam', help="...")
    parser.add_argument('--chrom', help="...", default=None)
    parser.add_argument('--start', help="...", type=int, default=None)
    parser.add_argument('--end', help="...", type=int, default=None)
    parser.add_argument('--out', help='...', default='out_phasogram.txt')
    return parser.parse_args(argv)


def main(argv):
    ''' docstring '''
    chrom = ''
    args = parse_args(argv)
    samfile = pysam.Samfile(args.bam, "rb")
    nextnuc = 0
    f_output = open(args.out, 'w')  # the output file
    buffer_dict = {True:  defaultdict(list), False:  defaultdict(list)}
    current_pos = {}
    for record in samfile.fetch(args.chrom, args.start, args.end):
        strand = record.is_reverse
        if record.mapq < _MINMAPQUALI:
            break  # do not analyze low quality reads
        elif record.tid != chrom:  # new chromosome
            chrom = record.tid
            buffer_dict[strand].clear()
            buffer_dict[not strand].clear()
            current_pos.clear()
            if strand:  # this is minus strand
                current_pos[strand] = record.aend  # because minus strand
                buffer_dict[strand]['start'].append(record.aend)
                buffer_dict[strand]['end'].append(record.pos)
                buffer_dict[strand]['score'].append(0)
            else:
                current_pos[strand] = record.pos  # plus strand
                buffer_dict[strand]['start'].append(record.pos)
                buffer_dict[strand]['end'].append(record.aend)
                buffer_dict[strand]['score'].append(0)

        else:
            if strand:  # this is minus strand
                current_pos[strand] = record.aend  # because minus strand
                buffer_dict[strand]['start'].append(record.aend)
                buffer_dict[strand]['end'].append(record.pos)
                buffer_dict[strand]['score'].append(0)
            else:
                current_pos[strand] = record.pos  # plus strand
                buffer_dict[strand]['start'].append(record.pos)
                buffer_dict[strand]['end'].append(record.aend)
                buffer_dict[strand]['score'].append(0)

            last_element = len(buffer_dict[strand]['start'])-1
            torefesh = []

            for idx in range(last_element-1):  # extra minus cause not check itself
                if buffer_dict[strand]['end'][last_element] >=  \
                        buffer_dict[strand]['start'][idx]:
                    buffer_dict[strand]['score'][last_element] += 1
                    buffer_dict[strand]['score'][idx] += 1
                if abs(buffer_dict[strand]['end'][idx] - current_pos[strand]) > _MAX_SIZE:
                    torefesh.append(idx)
            # print(len(torefesh),len(buffer_dict[strand]['start']),'laengde')
            # print(torefesh, 'before')

            for idx_out in torefesh:
                if buffer_dict[strand]['score'][idx_out] >= _MIN_DEPTH:
                    for i in range(last_element):
                        if not idx_out == i:
                            var = abs(buffer_dict[strand]['start'][idx_out] - buffer_dict[strand]['start'][i])
                            if var >= nextnuc:
                                print(var, file=f_output)
            nb_removed_already = 0

            for idx in torefesh:

                buffer_dict[strand]['start'].pop(idx - nb_removed_already)
                buffer_dict[strand]['end'].pop(idx - nb_removed_already)
                buffer_dict[strand]['score'].pop(idx - nb_removed_already)
                nb_removed_already += 1
    f_output.close()
    samfile.close()
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
