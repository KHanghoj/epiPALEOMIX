import pysam, re, sys, os, argparse
from pypeline.tools.commonutils import corr_chrom

def unpack(chrom, *rest):
    return chrom


def pa(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='Checking Bedfile chromosomes')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('BamPrefix', type=str)
    parser.add_argument('FastaPath', type=str)
    parser.add_argument('FastaPrefix', type=str)
    parser.add_argument('BedPath', type=str)
    return parser.parse_known_args(argv)

def get_fasta_chrom(args):
    faipath = args.FastaPath+'.fai'
    assert os.path.exists(faipath), ("reference fasta file must be indexed. "
                                     "Index with 'samtools faidx %s'" %
                                     args.FastaPath)
    with open(faipath, 'r') as fin:
        fastachroms = []
        for line in fin:
            fastachroms.append(unpack(*re.split(r'\s+', line.rstrip())))
    return fastachroms
    


def run(args):
    samfile = pysam.Samfile(args.BamPath, 'rb')
    bamchrom = [dic['SN'] for dic in samfile.header['SQ']]
    with open(args.BedPath, 'r') as fin:
        chromlast = ''
        bed_bam_chroms, bed_fasta_chroms = [], []
        for line in fin:
            chrom = unpack(*re.split(r'\s+', line.rstrip()))
            if chrom != chromlast:
                bed_bam_chroms.append(corr_chrom(args.BamPrefix, chrom))
                bed_fasta_chroms.append(corr_chrom(args.FastaPrefix, chrom))
                chromlast = chrom
        checked = [(c in bamchrom) for c in bed_bam_chroms]
        assert all(checked), \
            ('chromosome: "%s" in "%s" is/are not present in "%s"' %
             (', '.join(bed_bam_chroms[idx] for idx, c in enumerate(checked) if not c),
              os.path.basename(args.BedPath), os.path.basename(args.BamPath)))
        
        referencechrom = get_fasta_chrom(args)
        checked = [(c in referencechrom) for c in bed_fasta_chroms]
        assert all(checked), \
            ('chromosome: "%s" in "%s" is/are not present in "%s"' %
             (', '.join(bed_fasta_chroms[idx] for idx, c in enumerate(checked) if not c),
              os.path.basename(args.BedPath), os.path.basename(args.FastaPath)))


def main(argv):
    args, unknown = pa(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

