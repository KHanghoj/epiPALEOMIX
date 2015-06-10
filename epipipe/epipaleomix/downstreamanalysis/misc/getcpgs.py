from __future__ import print_function
#import os.path, re, pysam
import os.path, re
from epipaleomix.tools.commonutils import Cache
files = ['/home/krishang/data/bedfiles/methyl450k_1500_wochr.bed','/home/krishang/data/bedfiles/methyl450k_2000_wochr.bed']
referencepath='/home/krishang/data/reference_human/hs.build37.1.fa'
FMT = '{}\t{}\t{}\n'.format

def unpack(chrom, start, end, *rest):
    return str(chrom), int(start), int(end)


def getcpgs(seq, l):
    return (seq.count('C')+seq.count('G'))/float(l), seq.count('CG')


for fin in files:
    fasta = Cache(referencepath)
    with open(fin,'r') as infile:
        fout, _ = os.path.splitext(os.path.basename(fin))
        with open(fout+'.gccontentnew', 'w') as outfile:
            outfile.write(FMT('region','GCcontent', 'CpGcount'))
            for line in infile:
                chrom, start, end = unpack(*re.split(r'\s', line.rstrip()))
                regionname = '{}_{}_{}'.format(chrom,start, end)
                # start-1 because i wanna include the actual bases at position start.
                # bedfiles are always 1-based. fasta.fetch_string is 0based
                startzero = start-1
                outfile.write(FMT(regionname, *getcpgs(fasta.fetch_string(chrom, startzero, end-startzero), end-startzero)))



# def iterate_tree(tree, path=()):
#     if isinstance(tree, dict):
#         for key, value in tree.iteritems():
#             for elem in iterate_tree(value, path+(key,)):
#                 yield elem
#     else:
#         yield (path, tree)




# for fin in files:
#     curr, lst = '', []
#     with open(fin,'r') as infile:
#         fout, _ = os.path.splitext(fin)
#         with open(fout+'.gccontent', 'w') as outfile:
#             outfile.write(FMT('name','GCcontent', 'CpGcount'))
#             for line in infile:
#                 if line.startswith('>'):
#                     if curr:
#                         outfile.write(FMT('_'.join(re.split(r':|-',curr[1:])),
#                                           *getcpgs(''.join(lst))))
#                     curr = line.rstrip()
#                     lst = []
#                 else:
#                     lst.append(line.rstrip().upper())
