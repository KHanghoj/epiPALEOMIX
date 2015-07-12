from __future__ import print_function
import gzip, re, sys, argparse
from math import floor
FMT='{c}\t{s}\t{e}\t{depth}\t{score}\t{bedc}\n'.format

def unpackepipal(chrom, start, end, depth, score, bedc):
    return chrom, int(start), int(end), float(depth), float(score), str(bedc)

def unpack(chrom, start, end, depth, score, bedc):
    return {'c':str(chrom), 's':int(start), 'e':int(end),
            'depth':float(depth), 'score':float(score), 'bedc':str(bedc)}

def compareoutputfunc(fout, c, s, e, depth, score, bedc):
    mid=int(floor(e-s)+s)
    fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(c, mid-73, mid+73, depth, score, bedc))

def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk nucleosome data')
    parser.add_argument('infile', type=str)
    parser.add_argument('outfile', type=str)
    parser.add_argument('--mindyaddist', default=140, type=int)
    parser.add_argument('--minreaddepth', default=1, type=int)
    parser.add_argument('--scorecutoff', default=1, type=float)
    return parser.parse_args(args)

def nuclclean(args):
    MINDYADDIST = args.mindyaddist
    CUTOFF_READDEPTH = args.minreaddepth
    CUTOFF_SCORE = args.scorecutoff
    compareoutput = gzip.open('compareoutput.txt.gz', 'wb')
    with gzip.open(args.outfile, 'wb') as f_out:
        with gzip.open(args.infile, 'rb') as f_in:
            f_out.write(f_in.next())  # removes header
            lastdata = unpack(*re.split(r'\s+', next(f_in).rstrip()))
            
            while lastdata['score']<CUTOFF_SCORE or lastdata['depth']<CUTOFF_READDEPTH:
                lastdata = unpack(*re.split(r'\s+', next(f_in).rstrip()))

            for line in f_in:
                data = unpack(*re.split(r'\s+', line.rstrip()))
                if data['score']>=CUTOFF_SCORE and data['depth']>=CUTOFF_READDEPTH:                
                    if data['c'] != lastdata['c']:
                        f_out.write(FMT(**lastdata))
                        compareoutputfunc(compareoutput,**lastdata)
                        lastdata = data
                        continue
                    if data['s']-lastdata['e'] >= MINDYADDIST:
                        f_out.write(FMT(**lastdata))
                        compareoutputfunc(compareoutput,**lastdata)
                        lastdata = data
                    elif data['score']>lastdata['score']:
                        lastdata = data
                    # if this one overlap the previous or lest than dyaddist away. pick the best score. and assign data to lastdata
                    # if this one is more than dyaddist away or new chrom. print lastdata and assign
            f_out.write(FMT(**lastdata))
            compareoutputfunc(compareoutput,**lastdata)
    compareoutput.close()
def main(argv):
    args = p_a(argv)
    nuclclean(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


