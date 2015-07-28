import sys, re
import fileinput
from collections import namedtuple

OUTFMT = '{}\t{}\t{}\t{}\n'.format
BEDCOORD = '{}_{}_{}'.format
lasttup, intup = (), ()
chrom, start = '', 0
NA_TUP = namedtuple('row', 'c s e bedc')
# def parse_args(argv):
#     ''' docstring '''
#     parser = argparse.ArgumentParser()
#     parser.add_argument('inputfile', help='..', type=str)
#     parser.add_argument('mappauniqueness', help="..", type=float)
#     return parser.parse_known_args(argv)


## NA_TUP = namedtuple('row', 'c s e mappa bedc')

### bedtools intersect  -wb  -a ../reference_human/mappability_31/GENOME_31_50.40000-20000.mappability  -b conservednucleosomearray_HG19_wochr_withname.bed | sort -k 7,7n | python filtermappa_adv.py

# def unpack(c,s,e,midstuff,mappa,*rest):
#     return NA_TUP(c, int(s), int(e), BEDCOORD(*rest[:3]))
#     ## return NA_TUP(c, int(s), int(e), float(mappa), BEDCOORD(*rest[:3]))

def unpack(c,s,e,*rest):
    return NA_TUP(c, int(s), int(e), BEDCOORD(*rest[:3]))
    ## return NA_TUP(c, int(s), int(e), float(mappa), BEDCOORD(*rest[:3]))
    
for line in fileinput.input():
    intup = unpack(*re.split('\s+', line.rstrip('\n')))
    if not lasttup:
        start = intup.s
        chrom = intup.c
        lasttup = intup
    if intup.bedc == lasttup.bedc:
        if intup.s > lasttup.e:
            sys.stdout.write(OUTFMT(chrom, start, lasttup.e, intup.bedc))
            start = intup.s
    else:
        sys.stdout.write(OUTFMT(chrom, start, lasttup.e, lasttup.bedc))
        start = intup.s
        chrom = intup.c

    lasttup = intup

if intup:
    sys.stdout.write(OUTFMT(intup.c, start, intup.e, intup.bedc))

# if intup:
#     if start == intup.s and chrom==intup.c:  # if the very last line was on same bedc but nor overlapping
#         sys.stdout.write(OUTFMT(chrom, start, intup.e, intup.bedc))          
#     elif lasttup == intup:   ## the very last one was on 
#         sys.stdout.write(OUTFMT(intup.c, start, intup.e, intup.bedc))


# for line in fileinput.input():
#     intup = unpack(*line.rstrip('\n').split('\t'))
#     if intup.mappa >= mappauniqueness:
#         if not lasttup:
#             start = intup.s
#             chrom = intup.c
#             lasttup = intup
#         if intup.bedc == lasttup.bedc:
#             if intup.s > lasttup.e:
#                 sys.stdout.write(OUTFMT(chrom, start, lasttup.e, intup.bedc))
#                 start = intup.s
#         else:
#             sys.stdout.write(OUTFMT(chrom, start, lasttup.e, lasttup.bedc))
#             start = intup.s
#             chrom = intup.c
#         lasttup = intup

# if intup:
#     if start == intup.s and chrom==intup.c:  # if the very last is alone in same bedcoord
#         sys.stdout.write(OUTFMT(chrom, start, intup.e, intup.bedc))          
#     elif lasttup == intup:
#         sys.stdout.write(OUTFMT(intup.c, start, intup.e, intup.bedc))

        
# def parse_args(argv):
#     ''' docstring '''
#     parser = argparse.ArgumentParser()
#     parser.add_argument('inputfile', help='..', type=str)
#     parser.add_argument('mappauniqueness', help="..", type=float)
#     return parser.parse_known_args(argv)

            
# def main(argv):
#     args, unknown = parse_args(argv)
#     with open(args.inputfile, 'r') as infile:
#         lastbedcoord, lastline = '', []
#         chrom, start = '', ''
#         for line in infile:
#             inputline = line.rstrip('\n').split('\t')
#             if float(inputline[4]) >= args.mappauniqueness:
#                 # check if bedcoord is the same as before:
#                 bedcoord = BEDCOORD(*inputline[-3:])
#                 if bedcoord == lastbedcoord:
#                     lastline = inputline
#                     lastbedcoord = bedcoord
#                 else:
#                     end = lastline[2] if lastline[2]<lastline[7] else lastline[7]
#                     sys.stdout.write(OUTFMT(chrom, start, end, bedcoord))
#                     start = inputline[1] if inputline[1]> inputline[6] else inputline[6]
#                     chrom = inputline[0]
#                     lastline = inputline
#                     lastbedcoord = bedcoord
        

# if __name__ == '__main__':
#     sys.exit(main(sys.argv[1:]))
