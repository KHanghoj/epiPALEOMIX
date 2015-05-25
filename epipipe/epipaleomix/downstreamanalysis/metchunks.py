from __future__ import print_function
import gzip, re, sys, argparse
CHUNKRANGE = 2000
CUTOFF = 20
def unpackepipal(chrom_last, start_last, dea, tot, *pos):
    return chrom_last, int(start_last), int(dea), int(tot)

def unpackRRBS(chrom_last, start_last, *rest):
    return chrom_last, int(start_last), int(rest[-1])


def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk methylation data')
    parser.add_argument('metpath')
    parser.add_argument('--analysis', type=str, choices=['RRBS', 'epipaleomix'], default='epipaleomix')
    parser.add_argument('--cutoff', type=int, default=1)
    parser.add_argument('--chunksize', default=2000, type=int)
    return parser.parse_args(args)


def RRBS_anal(args):
    CHUNKRANGE = args.chunksize
    CUTOFF = args.cutoff
    with open(args.metpath) as f_in:
        checked = 0
        totdea, tottot, endval = 0, 0, CHUNKRANGE
        chrom_last = ''
        for line in f_in:
            if checked % 100000 == 0:
                print('checked {} methylated sites'.format(checked), file=sys.stderr)

            chrom, start, dea = unpackRRBS(*re.split(r'\s+', line.rstrip()))
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

def methylepipal(args):
    CHUNKRANGE = args.chunksize
    CUTOFF = args.cutoff
    with gzip.open(args.metpath) as f_in:
        print(f_in.next(), file=sys.stderr) # removes the header
        checked = 0
        totdea, tottot, endval = 0, 0, CHUNKRANGE
        chrom_last = ''
        for line in f_in:
            if checked % 100000 == 0:
                print('checked {} methylated sites'.format(checked), file=sys.stderr)

            chrom, start, dea, tot = unpackepipal(*re.split(r'\s+', line.rstrip()))
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

            if dea > 10:
                if (float(dea)/tot) == 1:
                    continue
            # may calc the means before summing to avoid a single CpG position skewing everything.
            # calculate the mean of means
            # tempmean += dea/tot
            totdea += dea
            tottot += tot
            checked += 1
        if tottot>CUTOFF:
            print(chrom_last, endval-CHUNKRANGE, endval, totdea, tottot,
                  float(totdea)/tottot, file=sys.stdout)


def main(argv):
    args = p_a(argv)
    run = methylepipal if args.analysis == 'epipaleomix' else RRBS_anal
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

