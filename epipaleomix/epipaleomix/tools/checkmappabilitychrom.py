import re, sys, os, argparse

def unpack(chrom, *rest):
    return chrom

def pa(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='Checking mappability chromosomes')
    parser.add_argument('Mappability', type=str)
    parser.add_argument('ChromUsed', type=list)
    return parser.parse_known_args(argv)

def checkasserts(checked, bed_bam_chroms, args):
    assert all(checked), \
        ('chromosome: "%s" in "%s" is/are not present in "%s"' %
         (', '.join(bed_bam_chroms[idx] for idx, c in enumerate(checked) if not c),
          os.path.basename(args.BedPath), os.path.basename(args.BamPath)))

def getchromprefix(c):
    """ removes number and uppercase characters from string to check prefix"""
    return ''.join(i for i in c if not (i.isdigit() or i.isupper()))

def checkmappabilityprefix(mappachrom, chromused):
    Mappaprefix = getchromprefix(mappachrom)
    Chrom_tobe_analyzed = getchromprefix(chromused)
    mappachrom = 'No prefix' if not Mappaprefix else Mappaprefix
    chromused = 'No prefix' if not Chrom_tobe_analyzed else Chrom_tobe_analyzed
    assert Mappaprefix == Chrom_tobe_analyzed, ("'--MappabilityPath' prefix: '%s' is "
                                                "different from 'ChromUsed'"
                                                " prefix '%s' used for GCcorrection") % (mappachrom, chromused)
    
def run(args):
    """ does not check if header is present in mappability file """ 
    with open(args.Mappability, 'r') as fin:
        chrom = unpack(*re.split(r'\s+', next(fin).rstrip()))
        for mappachrom in args.ChromUsed:
            checkmappabilityprefix(chrom, mappachrom)

def main(argv):
    args, unknown = pa(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

