import pysam, re, sys
from pypeline.tools.commonutils import corr_chrom

# bamfile input,
# bedfileinput
# output is true false
def unpack(chrom, *rest):
    return chrom

def pa():
    pass

#def main(argv):
def main():
#    args, unknown = pa(argv)
#    samfile = pysam.Samfile(args.bam, 'rb')
    samfile = pysam.Samfile('saqqaq_chrom22.bam', 'rb')
    bamchrom = [dic['SN'] for dic in samfile.header['SQ']]

    
#    with open(args.bed, 'r') as fin:
    with open('CTCF_human_2kb_rmoverlap.bed', 'r') as fin:
        chromlast = ''
        bedchroms = []
        for line in fin:
            chrom = unpack(*re.split(r'\s+', line.rstrip()))
            if chrom != chromlast:
#                bedchroms.append(corr_chrom(args.BamPrefix, chrom))
                bedchroms.append(corr_chrom('', chrom))
                chromlast = chrom
    checked = [(c in bamchrom) for c in bedchroms]
    assert all(checked), \
        ('chromosome: "%s" is/are not present in the bamfile' %
         ', '.join((bedchroms[idx] for idx, c in enumerate(checked) if not c)))


if __name__ == '__main__':
#    sys.exit(main(sys.argv[1:]))
    sys.exit(main())
