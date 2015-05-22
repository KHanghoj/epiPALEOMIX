from __future__ import print_function
import gzip, re, sys
CHUNKRANGE = 2000
CUTOFF = 10
def unpack(chrom_last, start_last, *rest):
    return chrom_last, int(start_last), int(rest[-1])


with open('ENCFF000LWB.bedMethyl.Bed') as f_in:
    checked = 0
    totdea, tottot, endval = 0, 0, CHUNKRANGE
    chrom_last = ''
    for line in f_in:
        if checked % 100000 == 0:
            print('checked {} methylated sites'.format(checked), file=sys.stderr)

        chrom, start, dea = unpack(*re.split(r'\s+', line.rstrip()))
        if chrom != chrom_last:
            if tottot>CUTOFF:
                print(chrom_last, endval-CHUNKRANGE, endval, totdea, tottot,
                      float(totdea)/tottot, file=sys.stdout)
            chrom_last = chrom
            totdea, tottot, endval = 0, 0, CHUNKRANGE
        if endval < start:
            if tottot>CUTOFF:
                print(chrom_last, endval-CHUNKRANGE, endval, totdea, tottot,
                      float(totdea)/tottot, file=sys.stdout)
                totdea, tottot = 0, 0
            while endval < start:
                endval += CHUNKRANGE
        totdea += dea
        tottot += 1
        checked += 1
    if tottot>CUTOFF:
        print(chrom_last, endval-CHUNKRANGE, endval, totdea, tottot,
              float(totdea)/tottot, file=sys.stdout)
