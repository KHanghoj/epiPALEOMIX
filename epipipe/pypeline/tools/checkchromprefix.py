import pysam, re
from pypeline.tools.commonutils import corr_chrom

# bamfile input,
# bedfileinput
# output is true false
def unpack(chrom, *rest):
    return chrom

def pa():
    pass

def main(argv):
    args, unknown = pa(argv)
    samfile = pysam.Samfile(args.bam, 'rb')
    bamchrom = [dic['SN'] for dic in samfile.header['SQ']]
    bedchroms = []
    chromlast = ''
    
    with open(args.bed, 'r') as fin:
        for line in fin:
            chrom = unpack(*re.split(r'\s+', line.rstrip()))
            if chrom != chromlast:
                bedchroms.append(corr_chrom(args.BamPrefix, chrom))
                chromlast = chrom
    checked = [(c in bamchrom) for c in bedchroms]
    assert all(checked), \
        ('chromosome: "%s" is/are not present in the bamfile' %
         ', '.join((bedchroms[idx] for idx, c in enumerate(checked) if not c)))
