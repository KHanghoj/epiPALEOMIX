from __future__ import print_function
import os.path, re
files = ['methyl450k_1500_wochr.fa', 'methyl450k_2000_wochr.fa']

FMT = '{}\t{}\t{}\n'.format


def getcpgs(seq):
    l = float(len(seq))
    return (seq.count('C')+seq.count('G'))/l, seq.count('CG')


for fin in files:
    curr, lst = '', []
    with open(fin,'r') as infile:
        fout, _ = os.path.splitext(fin)
        with open(fout+'.gccontent', 'w') as outfile:
            outfile.write(FMT('name','GCcontent', 'CpGcount'))
            for line in infile:
                if line.startswith('>'):
                    if curr:
                        outfile.write(FMT('_'.join(re.split(r':|-',curr[1:])),
                                          *getcpgs(''.join(lst))))
                    curr = line.rstrip()
                    lst = []
                else:
                    lst.append(line.rstrip().upper())
