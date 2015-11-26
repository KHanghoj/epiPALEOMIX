import pysam, argparse

def parse_args(argv):
    ''' docstring '''

    parser = argparse.ArgumentParser(prog='Get Min Max Readlengths')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('--MinMappingQuality', type=int)
    parser.add_argument('--MinAlignmentLength', type=int)
    parser.add_argument('--NoReadsChecked', type=int) 
    return parser.parse_known_args(argv)

def converttolist(argv):    
    arg = [argv.get('BamPath')]
    for key, pos in argv.items():
        if key.startswith('--'):
            arg.extend([key,str(pos)])
    return arg
    

## def main(pathtobam, MinMappingQuality, MinAlignmentLength, readswithinrange):
def main(argv):    
    args, unknown = parse_args(converttolist(argv))
    init = 1
    within_count = args.NoReadsChecked
    lst = []
    lstapp = lst.append
    with pysam.AlignmentFile(args.BamPath, "rb") as samfile:
        for rec in samfile:
            if (rec.alen < args.MinAlignmentLength or
                    rec.mapq < args.MinMappingQuality or
                    rec.is_unmapped):
                continue
            current=rec.alen
            lstapp(current)
            if init:
                upperbound = current
                lowerbound = current
                init = 0
                continue
            if current > upperbound:
                upperbound = current
                within_count = args.NoReadsChecked
            elif current < lowerbound:
                lowerbound = current
                within_count = args.NoReadsChecked
            else:
                within_count -= 1
            if not within_count:
                break
    lst.sort()
    top_ninetyfive = lst[int(0.95*len(lst))]
    return (lowerbound, upperbound, top_ninetyfive)




if __name__ == '__main__':
    pass
    ## print main('/Users/krishang/Desktop/example/saqqaq_chrom22.bam', 25, 27, 50000)
