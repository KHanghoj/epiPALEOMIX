from __future__ import print_function
import gzip, sys
from collections import namedtuple
NAtup = namedtuple('row', 'c, s d score relapos')
FMT='{}\t{}\t{}\n'.format

def splitbedc(bedc):
    ''' returns the start of bed region i.e. 1_10_20, return 10 '''
    return int(bedc.split('_')[1])

def unpack(c,s,d,score,bedc):
    s = int(s)
    return NAtup(c, s, float(d), float(score), s-splitbedc(bedc))

infile = sys.argv[1]
sizeofbedregion = 2020
with gzip.open(infile,'rb') as fin:
    next(fin)  ## removes header
    d = {k:[0,0] for k in xrange(0,sizeofbedregion+1)}
    for idx, line in enumerate(fin):
        if idx % 100000 == 0:
            sys.stderr.write(str(idx)+'\t')
        tup = unpack(*line.rstrip().split('\t'))
        curpos = d[tup.relapos]
        curpos[0] += tup.d
        curpos[1] += tup.score
sys.stderr.write('\n')
sys.stdout.write(FMT('relapos', 'depth', 'score'))
for k in sorted(d.iterkeys()):
    sys.stdout.write(FMT(k, d[k][0], d[k][1]))
