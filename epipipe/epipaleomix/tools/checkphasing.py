from __future__ import print_function
import gzip, re, sys
phaserange = 250

with gzip.open('gokcpgsnucleo.txt.gz') as f_in:
    totsum = 0
    counter = 0
    last_start = [0]
    for line in f_in:
        chrom, start, end, depth, score, bedcoord = re.split(r'\s+', line.rstrip())
        if float(score) > 1:
            if int(start)-last_start[-1] <= phaserange:
                if int(start)-last_start[-1] >= 140:
                    last_start.append(int(start))
                else:
                    pass # update the score if new score greater than last
            else:
                if len(last_start) > 2:
                    print(chrom, last_start[0], last_start[-1], len(last_start),
                          file=sys.stdout, sep='\t')
                    totsum += len(last_start)
                counter = 1
                last_start = [int(start)]
    if len(last_start) > 2:
        print(chrom, last_start[0], last_start[-1], len(last_start),
          file=sys.stdout, sep='\t')
#    print(totsum)
