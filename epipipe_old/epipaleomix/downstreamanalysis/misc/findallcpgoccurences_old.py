import pysam, re, sys
fastpath = 'chrome.fa'

CGpattern = re.compile('CG')
FMT = '{}\t{}\t{}\n'.format

def getindex(seq):
    for m in CGpattern.finditer(seq):
        yield m.start()

def getchroms(path):
    with open(path+'.fai') as f_chrom:
        return [line.split('\t', 1)[0] for line in f_chrom]
chroms = getchroms(fastpath)
for chrom in chroms:
    fasta = pysam.Fastafile(fastpath)
    size = fasta.get_reference_length(chrom)
    jumps = iter(xrange(1,size, int(1e7)))
    start = jumps.next()
    for end in jumps:
        for s in getindex(fasta.fetch('22',start,end)):
            beg = start+s+1  # goes form 0-based to 1-based to be used as BED file
            sys.stdout.write(FMT('22', beg, beg+1))
        start = end
sys.stdout.close()
