import pysam, re, sys, os, argparse
from pypeline.tools.commonutils import corr_chrom

# bamfile input,
# bedfileinput
# output is true false
def unpack(chrom, *rest):
    return chrom


def pa(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='Checking Bedfile chromosomes')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('bed', type=str)
    parser.add_argument('--FastaPrefix', dest='FastaPrefix', default='')
    parser.add_argument('--BamPrefix', dest='BamPrefix')
    return parser.parse_known_args(argv)


def main(argv):
    args, unknown = pa(argv)
    samfile = pysam.Samfile(args.BamPath, 'rb')
    bamchrom = [dic['SN'] for dic in samfile.header['SQ']]
    with open(args.bed, 'r') as fin:
        chromlast = ''
        bedchroms = []
        for line in fin:
            chrom = unpack(*re.split(r'\s+', line.rstrip()))
            if chrom != chromlast:
                bedchroms.append(corr_chrom(args.BamPrefix, chrom))
                chromlast = chrom
        checked = [(c in bamchrom) for c in bedchroms]
        assert all(checked), \
            ('chromosome: "%s" in "%s" is/are not present in "%s"' %
             (', '.join(bedchroms[idx] for idx, c in enumerate(checked) if not c),
              os.path.basename(args.bed), os.path.basename(args.BamPath)))

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
