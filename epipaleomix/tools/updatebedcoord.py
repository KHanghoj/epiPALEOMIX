import sys
import re
import fileinput
from collections import namedtuple

"""
takes a filter-by-uniqueness mappability file intersected with current bed file.
column 1-3 is the unique coordinates
column 4-6 is the original bed coordinates. They will still be used as bedcoordinates as they are unique for the input.
6+ is leftover from the original input.
"""
OUTFMT = '{}\t{}\t{}\t{}\n'.format
BEDCOORD = '{}_{}_{}_{strand}'.format


def unpack(c, s, e, *rest):
    strand = re.search(r"([+-])","".join(rest))
    if strand:
        bed = BEDCOORD(*rest[:3],
                       strand=strand.group(0))
    else:
        bed = BEDCOORD(*rest[:3],
                       strand='+')
    return NA_TUP(c, int(s), int(e), bed)
    # return NA_TUP(c, int(s), int(e), BEDCOORD(*rest[:3]))


lasttup, intup = (), ()
NA_TUP = namedtuple('row', 'c s e bedc')
chrom, start = '', 0


for line in fileinput.input():
    intup = unpack(*re.split('\s+', line.rstrip('\n')))
    if not lasttup:
        start = intup.s
        chrom = intup.c
        lasttup = intup
    if intup.bedc == lasttup.bedc:
        if intup.s > lasttup.e:
            sys.stdout.write(OUTFMT(chrom, start, lasttup.e, intup.bedc))
            start = intup.s
    else:
        sys.stdout.write(OUTFMT(chrom, start, lasttup.e, lasttup.bedc))
        start = intup.s
        chrom = intup.c
    lasttup = intup

if intup:
    sys.stdout.write(OUTFMT(intup.c, start, intup.e, intup.bedc))
else:
    raise RuntimeError("No Bed Regions overlap the uniquely mappable regions")
