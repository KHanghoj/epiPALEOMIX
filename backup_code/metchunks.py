from __future__ import print_function
import gzip, re, sys
CHUNKRANGE = 2000
CUTOFF = 20
def unpack(chrom_last, start_last, dea, tot, *pos):
    return chrom_last, int(start_last), int(dea), int(tot)


with gzip.open('BAM1_MethylMap_bed12345.txt.gz') as f_in:
    print(f_in.next(), file=sys.stderr) # removes the header
    checked = 0
    totdea, tottot, endval = 0, 0, CHUNKRANGE
    chrom_last = ''
    for line in f_in:
        if checked % 100000 == 0:
            print('checked {} methylated sites'.format(checked), file=sys.stderr)

        chrom, start, dea, tot = unpack(*re.split(r'\s+', line.rstrip()))
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

        # if (dea+1)/(tot+1) == 1 and dea > 9:
        #     print(chrom,start, dea, tot, 'likely a true SNV', file=sys.stderr)
        # else:
        #     totdea += dea
        #     tottot += tot
        if not ((dea+1)/(tot+1) == 1 and dea > 9):
            totdea += dea
            tottot += tot
        checked += 1
    if tottot>CUTOFF:
        print(chrom_last, endval-CHUNKRANGE, endval, totdea, tottot,
              float(totdea)/tottot, file=sys.stdout)


# with gzip.open('BAM1_MethylMap_bed12345.txt.gz') as f_in:
#     print(f_in.next(), file=sys.stderr) # removes the header
#     checked = 0
#     totdea, tottot = 0, 0
#     gen = iter(xrange(0, 9000000000, CHUNKRANGE))
#     start_last, end_last = gen.next(), CHUNKRANGE
#     chrom_last = ''
#     for line in f_in:
#         if checked % 100000 == 0:
#             print('checked {} methylated sites'.format(checked), file=sys.stderr)
#         chrom, start, dea, tot = unpack(*re.split(r'\s+', line.rstrip()))
#         if chrom != chrom_last:
#             if tottot:
#                 print(start_last, end_last, totdea, tottot,
#                       float(totdea)/(tottot+1), file=sys.stdout)
#             chrom_last = chrom
#             gen = iter(xrange(0, 9000000000, CHUNKRANGE))
#             start_last, end_last = gen.next(), CHUNKRANGE
#             totdea, tottot = 0, 0

#         if start_last+CHUNKRANGE < start:
#             if tottot:
#                 print(start_last, end_last, totdea, tottot,
#                       float(totdea)/(tottot+1), file=sys.stdout)
#                 totdea, tottot = 0, 0
#             while start_last+CHUNKRANGE < start:
#                 start_last = gen.next()
#                 end_last = start_last+CHUNKRANGE

#         if (dea+1)/(tot+1) == 1 and dea > 9:
#             print(chrom,start, dea, tot, 'likely a true SNV', file=sys.stderr)
#         else:
#             totdea += dea
#             tottot += tot
#         checked += 1

# with gzip.open('BAM1_MethylMap_bed12345.txt.gz') as f_in:
#     print(f_in.next(), file=sys.stderr) # removes the header
#     checked = 0
#     totdea, tottot = 0, 0
#     chrom_last, start_last, dea, tot = unpack(*re.split(r'\s+', f_in.next().rstrip()))
#     for line in f_in:
#         if checked % 100000 == 0:
#             print('checked {} methylated sites'.format(checked), file=sys.stderr)
#         chrom, start, dea, tot = unpack(*re.split(r'\s+', line.rstrip()))
#         if chrom == chrom_last and start_last + phaserange >= start:
#             if (dea+1)/(tot+1) == 1 and dea > 9:
#                 print(chrom,start, dea, tot, 'likely a true SNV', file=sys.stderr)
#             else:
#                 totdea += dea
#                 tottot += tot
#         else:
#             print(start_last, start_last+phaserange, totdea, tottot,
#                   float(totdea)/(tottot+1), file=sys.stdout)
#             totdea, tottot = 0, 0
#             chrom_last, start_last = chrom, start
#         checked += 1
