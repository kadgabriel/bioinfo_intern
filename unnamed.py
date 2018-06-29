import argparse, re, copy,  operator


def print_fragments(tempList):
	for i in range(0,len(tempList)):
		print("[%d] \"%s\"" % (i+1,tempList[i]))

def select_size(fragments, minsize, maxsize):
	selected_fragments = []
	for frag in fragments:
		if(len(frag) < maxsize and len(frag)>minsize):
			selected_fragments.append(frag)

	return selected_fragments

def digest(genome, p5, p3):
	fragments = re.split(p5+p3, genome)
	## Adding the hanging part
	curr_len = 0
	for i in range(0,len(fragments)-1):
		curr_len += len(fragments[i])

		temp_p5 = copy.copy(p5)
		temp_p3 = copy.copy(p3)

		## Add the 5' part
		if any(base in "[]" for base in temp_p5):
			temp_p5 = re.sub("\[.*?\]", genome[curr_len+temp_p5.find('[')], temp_p5)
			print(temp_p5)
		fragments[i] = fragments[i]+temp_p5
		curr_len += len(temp_p5)

		## Add the 3' part
		if any(base in "[]" for base in temp_p3):
			temp_p3 = re.sub("\[.*?\]", genome[curr_len+temp_p3.find('[')], temp_p3)
		fragments[i+1] = temp_p3+fragments[i+1]
	return fragments


def count_percent_unique(genome, fragmentList):
	count = 0
	for frag in fragmentList:
		count += len(frag)
	return count/float(len(genome))

def restriction_sites(enzyme):
	list_enzymes = {
		'SbfI'	:	'CCTGCA|GG',
		'ApeKI'	:	'G|CWGC',
		'EcoRI'	:	'G|AATTC',
		'PstI'	:	'CTGCA|G'
	}
	try:
		match_enzyme = list_enzymes[enzyme]
		if any(base in "NMRWYSKHBVD" for base in match_enzyme):
			match_enzyme = match_enzyme.replace("N", "[GCAT]")
			match_enzyme = match_enzyme.replace("M", "[CA]")
			match_enzyme = match_enzyme.replace("R", "[GA]")
			match_enzyme = match_enzyme.replace("W", "[AT]")
			match_enzyme = match_enzyme.replace("Y", "[CT]")
			match_enzyme = match_enzyme.replace("S", "[GC]")
			match_enzyme = match_enzyme.replace("K", "[GT]")
			match_enzyme = match_enzyme.replace("H", "[CAT]")
			match_enzyme = match_enzyme.replace("B", "[GCT]")
			match_enzyme = match_enzyme.replace("V", "[GCA]")
			match_enzyme = match_enzyme.replace("D", "[GAT]")
		site_p5, site_p3 = match_enzyme.split('|')
		return site_p5,site_p3
	except:
		print("nada")
		exit()

	exit()

def run_RE(enzyme):
	print()
	print(enzyme)
	p5,p3 = restriction_sites(enzyme)
	print(p5+p3)
	### READ FROM INPUT SEQUENCE FILE ###
	input_file  = open("input.txt", "r+")
	# genome ="GAGAGCTGCAGCGGCGCGGCAGCAA"
	for line in input_file:
		genome = line.strip()	
	# print(len(genome))
	## split genome according to RE (p5 and p3)
	fragments = digest(genome, p5, p3)

	## print the count of restriction sites
	print("Restriction sites:"+str(len(fragments)-1))

	# print all the fragments
	# print_fragments(fragments)

	# ## select the fragments based on size
	frag_select = select_size(fragments, minsize, maxsize)
	# # print_fragments(frag_select)
	print("Number of fragments filtered: ",end='')
	print(len(frag_select))

	# ## count selected fragments % in genome
	print("Percent coverage in genome: ",end='')
	print(count_percent_unique(genome,frag_select))
	return count_percent_unique(genome,frag_select)

def print_dict(items):
	for key in items:
   		print("Enzyme: %s\t Percent coverage: %f" % (key[0], key[1]))

if __name__ == '__main__':

	# restriction_sites('ApeKI')

	
	# print(p5+p3)
	minsize = int(input("Min fragment size: "))
	maxsize = int(input("Max fragment size: "))
	# minsize = 200
	# maxsize = 270
	list_enz = {}
	input_RE  = open("re.txt", "r+")
	for line in input_RE:
		enz = line.strip()
		list_enz[enz] = run_RE(enz)

	# print(list_enz)
	sorted_enz = sorted(list_enz.items(), key=operator.itemgetter(1),reverse=True)
	print_dict(sorted_enz)