import sys, re
import fileinput
from collections import namedtuple

OUTFMT = '{}\t{}\t{}\t{}\n'.format
BEDCOORD = '{}_{}_{}'.format
lasttup, intup = (), ()
chrom, start = '', 0
NA_TUP = namedtuple('row', 'c s e bedc')

def unpack(c,s,e,*rest):
    return NA_TUP(c, int(s), int(e), BEDCOORD(*rest[:3]))
    
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
