"""
RApyDS
Restriction Site Associated DNA Python-Digested Simulation

remove_repeat.py
"""

from __future__ import print_function
import sys, os
from operator import itemgetter

def remove_XAs(enzyme):
	if(os.path.isdir("bwa") == False and os.path.exists("bwa/aligned_pairs_"+enzyme+".sam")==False) :
		print("Run bwa_aln.sh first")
		raise SystemExit
	orig = open("bwa/aligned_pairs_"+enzyme+".sam","r")

	fragments = []
	XAs = []
	total = 0
	unique = 0
	repeat = 0
	for line in orig:
		if("LN:" in line):
			total = int(line.split("\t")[2][3:])
		if("Frag_" not in line):
			continue
		line = line.strip().rstrip().split("\t")
		frag_name = line[0].split("_")
		## all fragments's start and end are in this list fragments
		fragments.append([int(frag_name[2]),int(frag_name[3])])
		
		## get all XAs
		if(len(line) > 21):
			list_XA = line[21].split(":")[2].split(";")[:-1]
			for alt_hits in list_XA:
				alt = alt_hits.split(",")
				## each XA's location is in XAs list
				XAs.append(int(alt[1]))

		## get XT:A:U/R tags
		if(line[11] == 'XT:A:U'):
			unique += 1
		elif(line[11] == 'XT:A:R'):
			repeat +=1	

	orig.close()
	set_XAs = set(XAs) ## remove duplicates XAs
	fragments = sorted(fragments, key=itemgetter(0)) ## sort fragment's in increasing order
	ctr = 0
	for XA in set_XAs:
		if(XA < 0):
			XA += total
		for i in range(len(fragments)):
			if(fragments[i][0] > XA):
				break
			if(fragments[i][0] <= XA and XA <= fragments[i][1]):
				ctr+=1
				# print(XA)
				break

	#print(ctr)
	return ctr,unique,repeat

if __name__ == '__main__':
	pass
