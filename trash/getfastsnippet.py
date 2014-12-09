''' to get fasta.fa.fai: samtools faidx fasta.fa
This script creates a million bases from 10,000,000 to 11,000,000
of any given fasta file. in this case hs.build37 as the Labrana and saqqaq
is aligned against that build
'''


# import pysam, sys
# f = sys.argv[1]
# chrom = sys.argv[2]

# x = pysam.Fastafile(f)

# fasta_snip = open('fasta_snip.fa', 'w')
# y = x.fetch(chrom, start=10000000, end=11000000)
# fasta_snip.write(y)
# x.close()
# fasta_snip.close()




# def s(stril):
# 	del stril[:]

# sr = ['+']
# sr
# id(sr)
# s(sr)
# sr
# id(sr)




    # it = read_bed(args)
    # while True:  # while true needs to have a break within
                   # or if in a function, it needs a break
                   # outside the function if in a loop
                   # like below
    #     try:
    #         chrom, start, end = it.next()

    #         ....
    #     except StopIteration:
    #         break

    # def infinite():
    #     x = 0
    #     while True:
    #         yield x
    #         x += 1


string = 'CGTTTTCGASCG'
# string = 'TCGTCGT'
# pattern = 'CG'
# start  = 0
# while True:
#     idx = string.find(pattern,start)
#     if idx < 0:
#         break
#     print idx
#     start = idx + 1

def get_CG_index(string, pattern='CG'):
    start = 0
    while True:
        idx = string.find(pattern, start)
        if idx <0:
            break
        yield idx
        start = idx + 1
# string = 'TCsGTCsGT'
# a = get_CG_index(string)
# if list(get_CG_index(string, pattern='sdfsdf')) or list(get_CG_index(string, pattern='TCas')):
#     print True
# print list(a)
print(list(get_CG_index(string, pattern='CG')))
print(list(get_CG_index(string, pattern='TG')))
# print (a.next())
# print (a.next())
# print (a.next())
# import re
# ding = [m.start() for m in re.finditer('CG', string)]
# print (ding)
# ding = [m.start() for m in re.finditer('TC', string)]
# print (ding)
# # ding = [m for m in re.match('(?:CG|TC)', string)]
# # print (ding)

# ding = (re.finditer(r, string) for r in ['CG', 'TC'])
# # print ding.next()
# print(list(ding))


# for hit in get_CG_index(string):
#     print hit
