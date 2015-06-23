from __future__ import print_function
from scipy.stats import binom
import gzip, re, sys, argparse
FMT='{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format

MEANMETHYL = 0.07    # TpG/(CpG+(TpG-seqerror)) genome-wide for now
## THIS DOES NOT WORK MEANMETHYL = 0.4    # TpG/(CpG+(TpG-seqerror)) for sites containing a T only
SEQERROR = 0.003/3     # mean(c(708,495))*3/ (639816-(27566-mean(c(708,495))))
## divide with three as it is only the T error i'm interested in.

def unpackepipal(chrom, start, dea, tot, bedc):
    return chrom, int(start), int(dea), int(tot), str(bedc)


def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk methylation data')
    parser.add_argument('metpath')
    parser.add_argument('out')
    return parser.parse_args(args)


def writetofile(f_out, deaminatedsites, coverage, CpGsites, bedc):
    chrom, start, end = bedc.split('_')
    if deaminatedsites:
        methylprop = binom.cdf(deaminatedsites,coverage, MEANMETHYL)
    else:
        methylprop = binom.pmf(1,coverage, SEQERROR)*(1.0/coverage)  # the methylation proxy is the chance of having a sequencing error. observing a C when actual base was a T
    f_out.write(FMT(chrom, start, end, deaminatedsites, coverage, methylprop, CpGsites))


def methylepipal(args):
    with gzip.open(args.metpath, 'rb') as f_in:
        f_in.next() # removes the header
        checked = 0
        deaminatedsites, coverage, lastbedc = 0, 0, ''
        CpGsites = 0
        h = '#chrom\tpos\tend\tdeaminatedsites\tcoverage\tmethylprop\tCpGsites'
        with gzip.open(args.out, 'wb') as f_out:
            f_out.write(FMT(*h.split('\t')))
            for line in f_in:
                if checked % 100000 == 0:
                    print('checked {} methylated sites'.format(checked), file=sys.stderr)
                chrom, start, dea, tot, bedc = unpackepipal(*re.split(r'\s+', line.rstrip()))
                if bedc != lastbedc:
                    if coverage:
                        writetofile(f_out, deaminatedsites, coverage, CpGsites, lastbedc)
                    deaminatedsites, coverage, lastbedc  = 0, 0, bedc
                    CpGsites = 0
                if dea > 10:  # remove potential SNP's C>T. Does not represent methylation
                    if (float(dea)/tot) == 1:
                        continue
                deaminatedsites += dea
                coverage += tot
                CpGsites += 1
                checked += 1
            if coverage:
                writetofile(f_out, deaminatedsites, coverage, CpGsites, lastbedc)


def main(argv):
    args = p_a(argv)
    methylepipal(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

