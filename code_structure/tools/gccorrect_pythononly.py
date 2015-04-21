from __future__ import print_function
import re
import os
import sys
import numpy as np
# bamname = 'Saqqaq_GCcorrect'
# path = "/Users/kehanghoej/Desktop/GC"
# Desktop kehanghoej$ python gccorrect_pythononly.py GC Saqqaq


def getline(f):
    lst_read, lst_fasta = [], []
    with open(f, 'r') as fout:
        for line in fout:
            v, p, read, fasta = re.split(r'\s', line.rstrip())
            lst_read.append(float(read))
            lst_fasta.append(float(fasta))
        return v, np.asarray(lst_read), np.asarray(lst_fasta)+1


def calc(f):
    readlength, read, fasta = getline(f)
    grandmean = sum(read)/sum(fasta)
    dat = sum(abs(read / fasta - grandmean) * (fasta / sum(fasta)))
    return 0.5*dat/grandmean, readlength


def main(argv):
    path, bamname = argv
    valmax, optlength = 0, 0
    for dirpath, dirnames, fnames in os.walk(path):
        for fname in fnames:
            if fname.startswith(bamname):
                val, length = calc(os.path.join(dirpath, fname))
                if val > valmax:
                    valmax, optlength = val, length
    print(valmax, optlength, bamname, file=sys.stdout)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
