import argparse, re, copy,  operator

frag_location = []

def print_fragments(tempList):
	for i in range(0,len(tempList)):
		print("[%d] \"%s\"" % (i+1,tempList[i]))

def select_size(fragments, minsize, maxsize):
	## select the fragments given size
	selected_fragments = []
	for frag in fragments:
		if(len(frag) < maxsize and len(frag)>minsize):
			selected_fragments.append(frag)

	## create an array of fragment's location [0] start, [1] end
	global frag_location
	selected_locations = []
	for i in range(0,len(frag_location)):
		len_frag = frag_location[i][1] - frag_location[i][0]
		if(len_frag < maxsize and len_frag > minsize):
			selected_locations.append(frag_location[i])
	frag_location = selected_locations
	return selected_fragments

def digest(genome, p5, p3, p5_2, p3_2): 
	delimiters = [p5+p3, p5_2+p3_2]
	print (delimiters)
	regexPattern = '|'.join((delimiters))
	print (regexPattern)
	fragments = re.split(regexPattern, genome)

	## Adding the hanging part
	curr_len = 0
	global frag_location
	frag_location = []
	## Adding the hanging part
	for i in range(0,len(fragments)-1):
		temp_location = []
		temp_location.append(curr_len)
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
		
		temp_location.append(temp_location[0]+len(fragments[i]))
		frag_location.append(temp_location)

	temp_location = []
	len_genome = len(genome)
	temp_location.append(len_genome-len(fragments[len(fragments)-1]))
	temp_location.append(len_genome-1)
	frag_location.append(temp_location)
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


def parse_gff():
	input_gff = open("ecoli_sequence.gff3", "r+")

	gene_location = []
	for line in input_gff:
		parsed_line=line.strip().rstrip()
		if(len(parsed_line) > 0):
			parsed_line=parsed_line.split("\t")
			if(len(parsed_line) > 2):
				if(parsed_line[2] == "gene"):			## get only genes
					gene_location.append([int(parsed_line[3]),int(parsed_line[4])])	## add to gene_location the start and end
		
	return gene_location

def compare_gene(gene_location):
	gene_ctr = 0	## counter for gene location
	match_ctr = 0	## counter for number of matches (fragment in gene region)
	for i in range(0,len(frag_location)):
		## loop in the gene_location, fragment's location must be 
		## on the right of gene's start
		while(gene_location[gene_ctr][0] <= frag_location[i][0]):
			## return if gene_ctr exceeds
			if(gene_ctr == len(gene_location)-1):
				return match_ctr

			## if fragment is inside the gene region, we found a match!! else increment
			if(gene_location[gene_ctr][0] <= (frag_location[i][0]) and gene_location[gene_ctr][1] >= frag_location[i][1]):
				print(i)
				print("Fragment location\t"+str(frag_location[i]))
				print("Gene location\t\t"+str(gene_location[gene_ctr]))
				print()
				match_ctr += 1
				break
			else:
				gene_ctr += 1
				continue
	return match_ctr

def run_RE(enzyme, input_genome_file):
	print()
	print(enzyme)
	p5,p3 = restriction_sites(enzyme)
	print(p5+p3)
	### READ FROM INPUT SEQUENCE FILE ###
	input_file  = open(input_genome_file, "r+")

	genome=""
	for line in input_file:
		if(">" in line):	## skip lines with >
			continue
		genome += line.strip().rstrip()	## strip whitespaces
			
	## split genome according to RE (p5 and p3)
	global frag_location
	fragments = digest(genome, p5, p3,"","")

	## print the number of restriction sites
	print("Restriction sites:"+str(len(fragments)-1))

	# print all the fragments
	# print_fragments(fragments)

	# ## select the fragments based on size
	frag_select = select_size(fragments, minsize, maxsize)
	# # print_fragments(frag_select)
	print("Number of fragments filtered: ",end='')
	print(len(frag_select))

	## get all gene regions
	genes = parse_gff()

	## test case, if PstI, compare the fragments and genes
	if(enzyme == "PstI"):
		print(compare_gene(genes))
	# if(enzyme == "PstI"):
	# 	output = open("ecoli2.fasta", "w+")
	# 	for i in range(0,len(frag_select)):
	# 		output.write(">U00096.3 Fragment ")
	# 		output.write(str(i)+"\n")
	# 		output.write(frag_select[i])
	# 		output.write("\n")
	# 	output.close()

	## print all 
	# if(enzyme == "PstI"):
	# 	output = open("PstI_fragments", "w+")
	# 	for i in range(0,len(frag_select)):
	# 		output.write("Fragment ")
	# 		output.write(str(i)+"\n")
	# 		output.write(frag_select[i])
	# 		output.write("\n")
	# 	output.close()

	# ## count selected fragments % in genome
	print("Percent coverage in genome: ",end='')
	print(count_percent_unique(genome,frag_select))
	return count_percent_unique(genome,frag_select),len(frag_select)



def print_dict(items, items2):
	print("==============================")
	print("Enzyme\t Percent coverage\t Number of Fragments")
	for key in items:
   		print("%s\t %f\t\t %d" % (key[0], key[1], items2[key[0]]))


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='RADSeq python script')
	parser.add_argument('-i', nargs='?', default='ecoli.fasta')
	args = parser.parse_args()
	# print(args)
	# print(parser)

	# minsize = int(input("Min fragment size: "))
	# maxsize = int(input("Max fragment size: "))
	minsize = 500
	maxsize = 800
	percent_enz = {}
	num_enz={}
	input_RE  = open("re.txt", "r+")
	for line in input_RE:
		enz = line.strip()
		percent_enz[enz],num_enz[enz] = run_RE(enz, args.i)

	# print(list_enz)
	sorted_enz = sorted(percent_enz.items(), key=operator.itemgetter(1),reverse=True)
	print_dict(sorted_enz,num_enz)

	