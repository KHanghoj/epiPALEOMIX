from __future__ import print_function
import subprocess, sys, os.path
infile = sys.argv[1]
indata = subprocess.Popen(('zcat', infile), stdout=subprocess.PIPE)
next(indata.stdout,None)
lst = [int(line.rstrip().split('\t')[4]) for line in indata.stdout]
indata.kill()
outname = os.path.basename(infile).split('.')[0]
totaldatapoints = len(lst)
cov = 10
while len(lst)>(totaldatapoints/2):
   lst = [val for val in lst if val>=cov]
   cov+=10

sys.stdout.write('{}\t{}\n'.format(outname, cov))
