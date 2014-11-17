''' to get fasta.fa.fai: samtools faidx fasta.fa
This script creates a million bases from 10,000,000 to 11,000,000
of any given fasta file. in this case hs.build37 as the Labrana and saqqaq
is aligned against that build
'''


import pysam, sys
f = sys.argv[1]
chrom = sys.argv[2]

x = pysam.Fastafile(f)

fasta_snip = open('fasta_snip.fa', 'w')
y = x.fetch(chrom, start=10000000, end=11000000)
fasta_snip.write(y)
x.close()
fasta_snip.close()




def s(stril):
	del stril[:]

sr = ['+']
sr
id(sr)
s(sr)
sr
id(sr)
