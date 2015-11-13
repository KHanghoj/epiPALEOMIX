#!/usr/bin/env python

"""
Stripped down version of 'FilterUniqueBAM' by
:Author: Martin Kircher
:Contact: Martin.Kircher@eva.mpg.de
:Date: *08.10.2011
:Type: tool
:Input: BAM
:Output: BAM

Mark/Filter PCR duplicates for merged PE reads Reads BAM
from STDIN and writes BAM to STDOUT. All non-collapsed reads
as well as secondary/chinermic alignments, reads that have
failed QC and unmmaped reads, are written to STDOUT as is.

The input is assumed to be sorted by coordinates, and this
order is preservered, though individual reads at the same
position may be re-arranged).
"""

import sys
import pysam
from argparse import ArgumentParser


def calc_consensus(reads):
    count = len(reads)
    outread = None
    maxsumqual = 0
    for read in reads:
        nsum = sum(map(ord, read.qual))
        if nsum > maxsumqual:
            outread = read
            maxsumqual = nsum

        # LOOK FOR PREVIOUS PCR DUPLICATE COUNTS
        for key, value in read.tags:
            if key == "XP":
                count += value

    if not outread.tags:
        outread.tags = [("XP", count)]
    else:
        outread.tags = outread.tags + [("XP", count)]

    return outread


def get_consensus_se(reads):
    # DETERMINE MOST FREQUENT CIGAR LINE
    by_cigar = {}
    cigar_count = {}
    for read in reads:
        tcigar = tuple(read.cigar)
        if tcigar in by_cigar:
            cigar_count[tcigar] += 1
            by_cigar[tcigar].append(read)
        else:
            cigar_count[tcigar] = 1
            by_cigar[tcigar] = [read]

    to_sort = [(y, -len(str(x)), x) for (x, y) in cigar_count.iteritems()]
    to_sort.sort()
    selcigar = to_sort[-1][-1]
    reads = by_cigar[selcigar]

    return calc_consensus(reads)


def write_consensus_se(outfile, reads, remove_duplicates):
    consensus = get_consensus_se(reads)
    for read in reads:
        read.is_duplicate = (read is not consensus)
        if not (read.is_duplicate and remove_duplicates):
            outfile.write(read)


def _flush_buffer(outfile, curvariants, remove_duplicates):
    for value in curvariants.itervalues():
        write_consensus_se(outfile, value[0], remove_duplicates)
    curvariants.clear()


_FILTERED_FLAGS = 0x1     # PE reads
_FILTERED_FLAGS |= 0x4    # Unmapped
_FILTERED_FLAGS |= 0x100  # Secondary alignment
_FILTERED_FLAGS |= 0x200  # Failed QC
_FILTERED_FLAGS |= 0x800  # Chimeric alignment


def parse_args(argv):
    parser = ArgumentParser(usage="%(prog)s [options] < in.bam > out.bam")
    parser.add_argument("input", default="-", help="BAM file.", nargs="?")
    parser.add_argument("--remove-duplicates",
                        help="Remove duplicates from output; by default "
                             "duplicates are only flagged (flag = 0x400).",
                        default=False, action="store_true")

    return parser.parse_args(argv)


def main(argv):
    args = parse_args(argv)

    if args.input == "-" and sys.stdin.isatty():
        sys.stderr.write("STDIN is a terminal, terminating!\n")
        return 1
    elif sys.stdout.isatty():
        sys.stderr.write("STDOUT is a terminal, terminating!\n")
        return 1

    with pysam.Samfile(args.input, "rb") as infile:
        with pysam.Samfile("-", "wb", template=infile) as outfile:
            curpos = None
            curvariants = {}
            for (read_num, read) in enumerate(infile):
                if curpos and ((read.tid, read.pos) != curpos):
                    # Sort order is defined as ascending 'tid's and positions
                    if curpos > (read.tid, read.pos) and not read.is_unmapped:
                        sys.stderr.write("ERROR: Input file does not appear "
                                         "to be sorted by coordinates at "
                                         "record %i, aborting ...\n"
                                         % (read_num,))
                        return 1

                    _flush_buffer(outfile, curvariants,
                                  args.remove_duplicates)
                    curpos = None

                is_filtered = read.flag & _FILTERED_FLAGS
                is_collapsed = read.qname.startswith("M_")
                if is_filtered or not (read.qual and is_collapsed):
                    outfile.write(read)
                    continue

                curpos = (read.tid, read.pos)
                nkey = (read.is_reverse, read.pos, read.alen)
                if nkey in curvariants:
                    curvariants[nkey][0].append(read)
                    curvariants[nkey][1] += 1
                else:
                    curvariants[nkey] = [[read], 1]

            _flush_buffer(outfile, curvariants, args.remove_duplicates)

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
