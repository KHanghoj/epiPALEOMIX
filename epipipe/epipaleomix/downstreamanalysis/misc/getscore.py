from __future__ import print_function
import gzip, re, sys, argparse
from collections import deque
from itertools import islice
FMT='{c}\t{s}\t{score}\t{bedc}\n'.format
FMT='{}\t{}\t{}\t{}\n'.format
_SIZE = 147 
_OFFSET, _NEIGHBOR = 12, 25
_POSITION_OFFSET = _OFFSET+_NEIGHBOR
_TOTAL_WIN_LENGTH = _SIZE+(2*_OFFSET)+(2*_NEIGHBOR)
_HALFWINDOW =  (_TOTAL_WIN_LENGTH/2)
_CENTERINDEX = (_SIZE-1)/2
_SPACERMEAN = _NEIGHBOR+_NEIGHBOR

def unpack(c, s, depth, bedc):
    return {'c':str(c), 's':int(s), 'depth':float(depth), 'bedc':str(bedc)}

def calcscore(lst):
    spacer = (sum(islice(lst, 0, _NEIGHBOR)) + sum(islice(lst, _TOTAL_WIN_LENGTH-_NEIGHBOR, _TOTAL_WIN_LENGTH)))
    return lst[_NEIGHBOR+_OFFSET+_CENTERINDEX]-(spacer/float(_SPACERMEAN))

def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk nucleosome data')
    parser.add_argument('infile', type=str)
    parser.add_argument('outfile', type=str)
    return parser.parse_args(args)

def run(args):
    with gzip.open(args.outfile, 'wb') as f_out:
        with gzip.open(args.infile, 'rb') as f_in:
            #f_out.write(f_in.next())  # removes header
            maindepth = deque()
            for line in f_in:
                lastdata = unpack(*re.split(r'\s+', line.rstrip()))
                maindepth.append(lastdata['depth'])
                while len(maindepth)<_TOTAL_WIN_LENGTH:
                    lastdata = unpack(*re.split(r'\s+', next(f_in).rstrip()))
                    maindepth.append(lastdata['depth'])
                getscore = calcscore(maindepth)
                pos = lastdata['s']-_HALFWINDOW
                f_out.write(FMT(lastdata['c'], pos, getscore, lastdata['bedc']))
                maindepth.popleft()
                

def main(argv):
    args = p_a(argv)
    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


