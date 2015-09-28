from __future__ import print_function
import gzip, sys
from collections import namedtuple
NAtup = namedtuple('row', 'd score relapos bedc')
FMT='{}\t{}\t{}\n'.format

def splitbedc(bedc):
    ''' returns the start and end of bed region i.e. 1_10_20, return 10  20 '''
    c, s, e = bedc.split('_')
    return int(s), int(e)

def unpack(c,s,d,score,bedc):
    ''' do not need regend as in TSS as the orientation is not import for the CTCF '''
    regstart, regend = splitbedc(bedc)
    return NAtup(float(d), float(score), int(s)-regstart, bedc)

def calcbedrange(bedc):
    s,e = splitbedc(bedc)
    return int(e)-int(s)
    
# def unpack(c,s,d,score,bedc):
#     s = int(s)
#     return NAtup(c, s, float(d), float(score), s-splitbedc(bedc), bedc)


infile = sys.argv[1]
with gzip.open(infile,'rb') as fin:
    next(fin)  ## removes header
    init = True
    for idx, line in enumerate(fin):
        if idx % 500000 == 0:
            sys.stderr.write(str(idx)+" ")
        tup = unpack(*line.rstrip().split('\t'))
        if init:
            bedrange = calcbedrange(tup.bedc)  # need to do all this to initialize the range of the region as i do not want to hardcode
            d = {k:[0,0] for k in xrange(0,bedrange+1)}
            depthcorr = [0] * (bedrange+1)
            init = False
        curpos = d[tup.relapos]
        curpos[0] += tup.d
        curpos[1] += tup.score
        depthcorr[tup.relapos] += 1

sys.stderr.write('\n')
sys.stdout.write(FMT('relapos', 'depth', 'score'))
for k in sorted(d.iterkeys()):
    depth = d[k][0]/float(depthcorr[k])
    score = d[k][1]/float(depthcorr[k])
    sys.stdout.write(FMT(k, depth, score))

# sys.stderr.write('\n')
# sys.stdout.write(FMT('relapos', 'depth', 'score'))
# for k in sorted(d.iterkeys()):
#     sys.stdout.write(FMT(k, d[k][0], d[k][1]))
