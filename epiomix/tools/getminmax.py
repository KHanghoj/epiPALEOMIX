#!/usr/bin/env python
from __future__ import print_function
import pysam
import argparse


def parse_args(argv):
    ''' docstring '''

    parser = argparse.ArgumentParser(prog='Get Min Max Readlengths')
    parser.add_argument('BamPath', type=str)
    parser.add_argument('--MinMappingQuality', type=int, default=30)
    parser.add_argument('--MinAlignmentLength', type=int, default=25)
    parser.add_argument('--NoReadsChecked', type=int, default=5000) 
    return parser.parse_known_args(argv)

def converttolist(argv):    
    arg = [argv.get('BamPath')]
    for key, pos in argv.items():
        if key.startswith('--'):
            arg.extend([key,str(pos)])
    return arg
    

def main(argv):    
    args, unknown = parse_args(converttolist(argv))
    # args, unknown = parse_args(argv)  # for testing
    lst = []
    noreads = args.NoReadsChecked
    lstapp = lst.append
    with pysam.AlignmentFile(args.BamPath, "rb") as samfile:
        while noreads:
            try:
                rec = samfile.next()
            except StopIteration:
                break
            if (rec.alen < args.MinAlignmentLength or
                rec.mapq < args.MinMappingQuality or
                rec.is_duplicate or
                rec.is_secondary or      # this is primarily for BWA MEM
                rec.is_supplementary or  # this is primarily for BWA MEM
                rec.is_qcfail or
                rec.is_unmapped):
                continue
            noreads -= 1
            lstapp(rec.alen)
    lst.sort()
    top_ninetyfive = lst[int(0.95*len(lst))]
    return (min(lst), max(lst), top_ninetyfive)


if __name__ == "__main__":
    import sys
    sys.exit(main(sys.argv[1:]))
