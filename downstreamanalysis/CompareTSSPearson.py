# from __future__ import print_function
# import gzip, sys, re
# from itertools import izip
# from scipy import signal
# from copy import copy
# from collections import namedtuple
# NAtup = namedtuple('row', 'd score relastart bedc')
# FMT='{}\t{}\t{}\t{}\n'.format

# def splitbedc(bedc):
#     ''' returns the start and end of bed region i.e. 1_10_20, return 10  20 '''
#     c, s, e = bedc.split('_')
#     return int(s), int(e)

# def unpack(c,s,d,score,bedc):
#     s = int(s)
#     regstart, regend = splitbedc(bedc) 
#     return NAtup(float(d), float(score), s-regstart, bedc)

# def getstrand(f):
#     dic = {}
#     with open(f,'r') as fbed:
#         for line in fbed:
#             pos, strand, rest = line.rstrip().split('\t', 2)
#             dic[pos] = strand
#     s,e = splitbedc(pos)
#     return dic, e-s

# def retrievelst(lst, strand):
#     beg, end = CENTERDIC[strand]
#     return calcfourier(lst[beg:end])

# def calcfourier(lst):
#     freq, spec = signal.periodogram(lst)
#     cov = sum(1 for x in lst if x!=0)/float(len(lst))
#     return (cov, 1/freq[1:], spec[1:])

# def writetofile(lastbed, cov, freq, spec):
#     sys.stdout.write(FMT(lastbed, cov, freq[spec.argmax()], spec.max()))
#     # for f, s in izip(freq, spec):
#     #     sys.stdout.write(FMT(lastbed, cov, f, s))

# infile = sys.argv[1]
# # if re.search('TSS', infile):
# #     infolist = '/home/krishang/data/bedfiles/TSS_join_to_epipaleomix.txt'
# # elif re.search('PROM', infile):
# #     infolist = '/home/krishang/data/bedfiles/PROM_join_to_epipaleomix.txt'
# # else:
# #     raise IOError, 'Input file is unknown; takes either TSS or PROM file'
# if re.search('2kb', infile):
#     infolist = "/home/krishang/data/bedfiles/TSS_join_to_epipaleomix.txt"
# elif re.search('15', infile):
#     infolist = "/home/krishang/data/bedfiles/PROM_join_to_epipaleomix.txt"
# else:
#     raise IOError, "Input file is unknown; takes either TSS or PROM file"

# stranddic, bedrange = getstrand(infolist)
# lsttemplate = [0] * (bedrange+1)
# if bedrange == 2000:
#     CENTERDIC  = {'+':[800, 1800], '-':[200,1200]}
# elif bedrange == 1400:
#     CENTERDIC  = {'+':[200, 1200], '-':[200,1200]}

# with gzip.open(infile,'rb') as fin:
#     next(fin)  ## removes header
#     lastbed = ""
#     d = []
#     sys.stdout.write(FMT('region', 'cov', 'peak', 'value'))
#     for idx, line in enumerate(fin):
#         if idx % 500000 == 0:
#             sys.stderr.write(str(idx)+" ")
#         tup = unpack(*line.rstrip().split("\t"))
#         if tup.bedc != lastbed:
#             curr_strand = stranddic[tup.bedc]
#             if d:
#                 writetofile(lastbed, *retrievelst(d, stranddic[lastbed]))
#             d = copy(lsttemplate)  # need to make new object every time
#         d[tup.relastart] += tup.score
#         lastbed = tup.bedc

# sys.stderr.write('\n')
# writetofile(lastbed, *retrievelst(d, stranddic[lastbed]))

