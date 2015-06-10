from __future__ import print_function
import sys, os.path, re
#files = ['PROM_minusstrand_autosom_wochr.fa', 'PROM_plusstrand_autosom_wochr.fa']
files = ['PROM_autosom_wochr.fa']
WINDOWLENGTH = 500
HIGHCPG = 0.75
MIDCPG = 0.48
HIGHCGs = 0.55
FMT = '{}\t{}\n'.format


def CpGyield(seq):
    l = len(seq)
    for idx in xrange(l-WINDOWLENGTH+1):
        dat = seq[idx:(idx+WINDOWLENGTH)]
        cgs = (dat.count('C')+dat.count('G'))/float(WINDOWLENGTH)
        exp = (cgs/2)*(cgs/2)
        obs = dat.count('CG')/float(WINDOWLENGTH)
        yield cgs, obs/exp

def getcpgs(seq):
    flag = 'LOW'
    for cgs, cpg in CpGyield(seq):
        if cgs >= HIGHCGs and cpg >= HIGHCPG:
            return 'HIGH', cgs, cpg
        elif cpg >= MIDCPG:
            flag = 'INTERMEDIATE'
    return flag

for fin in files:
    curr, lst = '', []
    with open(fin,'r') as infile:
        fout, _ = os.path.splitext(fin)
        with open(fout+'.promcat', 'w') as outfile:
            outfile.write(FMT('name','PromGroup'))
            for line in infile:
                if line.startswith('>'):
                    if curr:
                        outfile.write(FMT('_'.join(re.split(r':|-',curr[1:])),
                                          getcpgs(''.join(lst))))
                    curr = line.rstrip()
                    lst = []
                else:
                    lst.append(line.rstrip().upper())
