from __future__ import print_function
import gzip, re, sys, argparse
FMT='{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format

def unpackepipal(chrom, start, dea, tot, bedc):
    return chrom, int(start), int(dea), int(tot), str(bedc)


def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk methylation data')
    parser.add_argument('metpath')
    parser.add_argument('out')
    return parser.parse_args(args)


def writetofile(f_out, totdea,tottot, totmean, meanchecked, bedc):
    chrom, start, end = bedc.split('_')
    f_out.write(FMT(chrom, start, end, totdea, tottot,
                         float(totdea)/tottot, totmean, meanchecked, float(totmean)/meanchecked))


def methylepipal(args):
    with gzip.open(args.metpath, 'rb') as f_in:
        f_in.next() # removes the header
        checked = 0
        totdea, tottot, lastbedc = 0, 0, ''
        totmean, meanchecked = 0, 0
        h = '#chrom\tpos\tend\tdea\ttotalcoverage\tdearate\tsumofmeans\tCpGsites\tmeanofmeans'
        with gzip.open(args.out, 'wb') as f_out:
            f_out.write(FMT(*h.split('\t')))
            for line in f_in:
                if checked % 100000 == 0:
                    print('checked {} methylated sites'.format(checked), file=sys.stderr)
                chrom, start, dea, tot, bedc = unpackepipal(*re.split(r'\s+', line.rstrip()))
                if bedc != lastbedc:
                    if tottot:
                        writetofile(f_out, totdea, tottot, totmean, meanchecked, lastbedc)
                    totdea, tottot, lastbedc  = 0, 0, bedc
                    totmean, meanchecked = 0, 0
                if dea > 10:  # remove potential SNP's C>T. Does not represent methylation
                    if (float(dea)/tot) == 1:
                        continue
                totdea += dea
                tottot += tot
                totmean += float(dea)/tot
                meanchecked += 1
                checked += 1
            if tottot:
                writetofile(f_out, totdea, tottot, totmean, meanchecked, lastbedc)


def main(argv):
    args = p_a(argv)
    methylepipal(args)
#    run = methylepipal if args.analysis == 'epipaleomix' else RRBS_anal
#    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

