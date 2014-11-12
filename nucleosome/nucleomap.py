#!/opt/local/bin/python
from __future__ import print_function

#### imports
import numpy as np
#from fileinput import inpua
from sys import exit, argv
from os import path
import pysam
import collections


#### Constants:
size = 147 # the window
offset = 12 # offset between spacer and nucleosomal DNA
neighbor = 25 #size of flanking regions to be considered for the log-odd ration score 

sites = 0
windows = [] # the sum of 25+12+147+12+25=221 
positions = [] # the chromosomal position of the nuclesome
basecomposi = [] # single nucleotide base composition
posi=[]
cnt=collections.Counter() # this is for counting bases in basecomposi

### unixcommand:
#$python mymapnucleosome1_pysam1.py test.bam




''' 
Object: call the nucleosome center if exists and return key=POSITION, value=DEPTH
It extends the values around the centre of the nucleosome if identical maximal depth
'''
a= sadf 
asdf= sad
def call_max(sliceT):
	call = {} # this is the return dic with position and depth
	maxdepth = max(sliceT)
	if maxdepth <=2: # do not take maxseqdep of 2 or lower into account
		call['NA'] = 1
		return call
	center_index = (len(sliceT)-1)/2
	# print(maxdepth, 'maxdepth', sliceT[center_index],'center',center_index)
	if sliceT[center_index] == maxdepth:
		call[center_index] = sliceT[center_index]
		i=1
		while (sliceT[(center_index - i)] == maxdepth) & (i < center_index):
			call[(center_index-i)] = sliceT[(center_index-i)] ## extends 5'
			i+=1
			# print(i,'I')
		i=1
		while (sliceT[(center_index + i)] == maxdepth) & (i < center_index):
			call[(center_index+i)] = sliceT[(center_index+i)] ## extends 3'
			i+=1
	else: call['NA']=1 # not callable. i.e the centre is not maximal
    return call 


def call_flanking(sliceT):
	maxdepth = max(sliceT)
	center_index = (len(sliceT)-1)/2
	if sliceT[center_index] == maxdepth:
		return sliceT[center_index]
	else:
		return 1

### here starts the script:

total_win_len = size+(2*offset)+(2*neighbor)

if not path.exists(argv[1]):
    exit('ERROR: Database %s was not found!' % argv[1])

f = argv[1]
samfile = pysam.Samfile( f, "rb" )
# result = list(x.n for x in samfile.pileup())


for pileupcolumn in samfile.pileup('22',16050360,16057100):
	s_depth = pileupcolumn.n # get the depth
	posi = list((pileupcolumn.tid+1,pileupcolumn.pos+1)) # creates a tuple of chr,pos
	
	sites+=1

	base = [x.alignment.seq[x.qpos] for x in pileupcolumn.pileups][0] # get the consensus nucleotide at each site
	if (sites<= total_win_len):
		windows.append(int(s_depth)) # the actual windows with pileup
		positions.append(posi) # positions
		

		# all commented is for base composition
		# basecomposi.append(base)

	else:
		windows.pop(0)
		windows.append(int(s_depth))
		positions.pop(0)
		positions.append(posi)
		# basecomposi.pop(0)
		# basecomposi.append(base)

	if len(windows) == total_win_len:
		calls_center = call_max(windows[neighbor+offset:neighbor+
			offset+size])

		if not 'NA' in calls_center.keys(): # checks if the center of nuclesome == maxdepth
			calls_spacerL = call_flanking(windows[0:neighbor])
			calls_spacerR = call_flanking(windows[neighbor-1+offset+size+offset:total_win_len-1])
			
			mean_spacer = 1.0 + 0.5 * (calls_spacerL + calls_spacerR)

			# cnt.update(''.join(basecomposi)) # concates the basecompi list into a string
			# tot = sum(cnt.values()) 
			# print([[k,(float(v)/tot)] for k,v in cnt.items()])
			# print(cnt)
			# cnt.clear()
			
			for key, value, in calls_center.items():
				score = np.log(float(value)/mean_spacer)
				position_call = positions[(neighbor+offset+int(key))]

				# print(position_call[0],position_call[1],value,score)
				print('{0[0]},{0[1]},{1},{2}'.format(position_call,value,score))

	

samfile.close()


a = []

# from collections import defaultdict
# def C(s):
# 	dinucleot = 2
# 	d = defaultdict(int)
# 	# d = {}
# 	for i in xrange(len(s)-dinucleot-1):
# 		d[s[i:i+dinucleot]] += 1
# 	return d

# print(C('ACGTAGCTACGACTCGATGCATGTAGCTAGCTAGCTACGTCAGTAGCTGACTGATCGATCGTCGTAGCTGACTGATCGATCGATCGAATTCCGG'))
