from __future__ import print_function
import re
## zcat ipi.HUMAN.v3.01.dat.gz | grep '//\|^AC\|REFSEQ'> tmp.file

s = open('IPIgenelist.txt','r')
l = set([l.rstrip() for l in s ])
s.close()
FMT='{}\t{}\n'.format
with open('output.txt','w') as fout:
    with open('tmp.file','r') as fin:
        for line in fin:
            if line.startswith('AC'):
                ipi = re.split(r'\s+',line.rstrip('\n|;'))[1]
                if ipi in l:
                    try:
                        # the (\d+) is used if group 1 is needed, i.e. the numbersonly. remember the r as we have a \ in the expression
                        npname = re.search(r'[NX][P]_(\d+)',next(fin)).group(0)
                        fout.write(FMT(ipi,npname))
                    except AttributeError:
                        pass

