from __future__ import print_function
import gzip, sys, re
from collections import namedtuple
NAtup = namedtuple('row', 'd score relastart relaend bedc')
FMT='{}\t{}\t{}\n'.format

def splitbedc(bedc):
    ''' returns the start and end of bed region i.e. 1_10_20, return 10  20 '''
    c, s, e = bedc.split('_')
    return int(s), int(e)

def unpack(c,s,d,score,bedc):
    s = int(s)
    regstart, regend = splitbedc(bedc) 
    return NAtup(float(d), float(score), s-regstart,  regend-s, bedc)

def getstrand(f):
    dic = {}
    with open(f,'r') as fbed:
        for line in fbed:
            pos, strand, rest = line.rstrip().split('\t', 2)
            dic[pos] = strand
    return dic

stranddic = getstrand('/home/krishang/data/bedfiles/TSS_join_to_epipaleomix.txt')
infile = sys.argv[1]
sizeofbedregion = 2000
with gzip.open(infile,'rb') as fin:
    next(fin)  ## removes header
    d = {k:[0,0] for k in xrange(0,sizeofbedregion+1)}
    lastbed = ""
    for idx, line in enumerate(fin):
        if idx % 500000 == 0:
            sys.stderr.write(str(idx)+" ")
        tup = unpack(*line.rstrip().split("\t"))
        if tup.bedc != lastbed:
            curr_strand = stranddic[tup.bedc]
        if curr_strand == '+':
            curpos = d[tup.relastart]
        else:
            curpos = d[tup.relaend]
        curpos[0] += tup.d
        curpos[1] += tup.score
        lastbed = tup.bedc
sys.stderr.write('\n')
sys.stdout.write(FMT('relapos', 'depth', 'score'))
for k in sorted(d.iterkeys()):
    sys.stdout.write(FMT(k, d[k][0], d[k][1]))
