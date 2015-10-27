import re, sys, os, argparse

def unpack(chrom, *rest):
    return str(chrom)

def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='Checking mappability chromosomes')
    parser.add_argument('Mappability', type=str)
    parser.add_argument('ChromUsed', type=list)
    return parser.parse_known_args(argv)

def getmappachroms(args):
    with open(args.Mappability, 'r') as fin:
        mappa = set()
        lastc = ''
        for line in fin:
            mappa_chrom = re.split(r'\s+', line.rstrip())[0]
            if mappa_chrom != lastc:
                mappa.add(mappa_chrom)
    return mappa


def run(args):
    """ does not check if header is present in mappability file """ 
    mappachroms = getmappachroms(args)
    print(args.ChromUsed)
    for chromused in args.ChromUsed:
        assert (chromused in mappachroms) , ("'ChromUsed' chromosome: '%s' chosen for GC-correction is not"
                                           " available in the Mappability file (--MappabilityPath):"
                                           "\n\t'%s'\nMake sure to just same prefixes"
                                           "in ChromUsed and Mappability file\n"
                                           "ChromUsed prefixes can be changed in the makefile easily") % (chromused, ' '.join(sorted(mappachroms)))

def main(argv):
    args, unknown = parse_args(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

