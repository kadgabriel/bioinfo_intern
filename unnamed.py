#!/usr/bin/python3
from __future__ import print_function
import argparse, re, copy, operator
import sys, time

def print_fragments(tempList):
	for i in range(0,len(tempList)):
		print("[%d] \"%s\"" % (i+1,tempList[i]))

def print_dict(items, items2):
	print("==============================")
	print("Enzyme\t Percent coverage\t Number of Fragments")
	for key in items:
   		print("%s\t %f\t\t %d" % (key[0], key[1], items2[key[0]]))

##	FUNCTION TO DIGEST A GENOME GIVEN P5 AND P3 SITES ##
def digest(genome, p5, p3):
	fragments = re.split(p5+p3, genome)		## split the genome into fragments
	
	curr_len = 0				## temporary holder for current base in genome
	new_fragments = []			## temporary list holder for digested fragments
	## Adding the hanging part	
	for i in range(0,len(fragments)-1):
		temp_frag = []			## temporary list holder for current fragment
		temp_start = curr_len	## take not of curr_len, it will be the start location of the fragment
		curr_len += len(fragments[i])

		temp_p5 = copy.copy(p5)
		temp_p3 = copy.copy(p3)

		## Add the 5' part
		if any(base in "[]" for base in temp_p5):
			temp_p5 = re.sub("\[.*?\]", genome[curr_len+temp_p5.find('[')], temp_p5)
			# print(temp_p5)
		fragments[i] = fragments[i]+temp_p5
		curr_len += len(temp_p5)

		## Add the 3' part
		if any(base in "[]" for base in temp_p3):
			temp_p3 = re.sub("\[.*?\]", genome[curr_len+temp_p3.find('[')], temp_p3)
		fragments[i+1] = temp_p3+fragments[i+1]
		temp_frag.append(fragments[i])
		temp_frag.append(temp_start)
		temp_frag.append(temp_start+len(fragments[i]))
		new_fragments.append(temp_frag)

	## append last fragment to list
	len_genome = len(genome)
	temp_frag = []
	temp_frag.append(fragments[len(fragments)-1])
	temp_frag.append(len_genome-len(fragments[len(fragments)-1]))
	temp_frag.append(len_genome-1)
	new_fragments.append(temp_frag)
	return new_fragments

##	FUNCTION TO SELECT ONLY FRAGMENTS WITHIN A GIVEN SIZE
def select_size(fragments, minsize, maxsize):
	## select the fragments given size
	selected_fragments = []
	for frag in fragments:
		## fragment's start and end must be inside the gene region
		if(len(frag[0]) < maxsize and len(frag[0]) > minsize):
			selected_fragments.append(frag)
	return selected_fragments

##	FUNCTION TO PARSE THE ENZYME DATABASE 
##	format:	each enzyme in separate lines
##		ex.	SbfI,CCTGCA|GG
##			ApeKI,G|CWGC	
def parse_enzymedb(enzyme_db_file):
	list_enzymes = {}
	input_db = open(enzyme_db_file, "r+")
	line_no = 1
	for line in input_db:
		line = line.strip().rstrip().split(",")		## strip strip and split

		## catch 'em all errors
		if(len(line)!=2):
			print("Error in restriction enzyme database file line no "+str(line_no))
			raise SystemExit
		else:
			search=re.compile(r'[GCATNMRWYSKHBVD]+[|]+[GCATNMRWYSKHBVD]+').search
			if(bool(search(line[1])) == False):
				print("Error in restriction enzyme database file line no "+str(line_no)+". Invalid letters in sequence.")
				raise SystemExit
			else:
				line_no+=1
				list_enzymes[line[0]] = line[1]

	if(len(list_enzymes) < 0):
		print("No restriction enzymes found in "+enzyme_db_file)
		raise SystemExit

	return list_enzymes


##	PARSER TO ENZYME
def restriction_sites(enzyme, list_enzymes):

	try:	## test if the RE is in loaded DB
		match_enzyme = list_enzymes[enzyme]
	except:
		print("Restriction enzyme "+enzyme+" is not in database")
		raise SystemExit

	## replace wildcards
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
	
	## split into p5 and p3
	site_p5, site_p3 = match_enzyme.split('|')
	return site_p5,site_p3


##	FUNCTION TO PARSE THE GENE ANNOTATION
##	format:	GFF
def parse_gff(annotation_file):
	input_gff = open(annotation_file, "r+")

	gene_location = []
	for line in input_gff:
		parsed_line=line.strip().rstrip()	## strip whitespaces
		if(len(parsed_line) > 0):
			parsed_line=parsed_line.split("\t")
			if(len(parsed_line) > 2 and parsed_line[0] != '#'):
				## get only genes for now
				if(parsed_line[2] == "gene"):			
					gene_location.append([int(parsed_line[3]),int(parsed_line[4])])	## add to gene_location the start and end
	## exit if there are no genes parsed
	if(len(gene_location) == 0):
		print("Check the annotation file. Must be in GFF format or contains none of the features wanted.")
		raise SystemExit
	return gene_location


##	Function that counts how many fragments are within the gene region
def compare_gene(gene_location, fragments):
	gene_ctr = 0	## counter for gene location
	match_ctr = 0	## counter for number of matches (fragment in gene region)
	for i in range(0,len(fragments)):
		## loop in the gene_location, fragment's location must be 
		## on the right of gene's start
		while(gene_location[gene_ctr][0] <= fragments[i][1]):
			## return if gene_ctr exceeds
			if(gene_ctr == len(gene_location)-1):
				return match_ctr

			## if fragment is inside the gene region, we found a match!! else increment
			if(gene_location[gene_ctr][0] <= (fragments[i][1]) and gene_location[gene_ctr][1] >= fragments[i][2]):
				# print(i)
				# print("Fragment location\t"+str([fragments[i][1],fragments[i][2]]))
				# print("Gene location\t\t"+str(gene_location[gene_ctr]))
				# print()
				match_ctr += 1
				break
			else:
				gene_ctr += 1
				continue
	return match_ctr


def run_RE(enzyme, parsed, args):

	print()
	# print(enzyme)
	print(enzyme,end='\t')
	
	p5,p3 = restriction_sites(enzyme,parsed['db'])
	# print(p5+p3)
	# print(p5+p3,end='\t')

	### READ FROM INPUT SEQUENCE FILE ###
	input_file  = open(args.i, "r+")

	genome=""
	for line in input_file:
		if(">" in line):	## skip lines with >
			continue
		genome += line.strip().rstrip()	## strip whitespaces
			
	## split genome according to RE (p5 and p3)
	fragments = digest(genome, p5, p3)

	## print the number of restriction sites
	print(str(len(fragments)-1),end='\t')
	# print("Restriction sites:"+str(len(fragments)-1))

	# ## select the fragments based on size
	frag_select = select_size(fragments, args.min, args.max)

	# print("Number of fragments filtered: ",end='')
	# print(len(frag_select))
	print(len(frag_select),end='\t')

	## get all gene regions and 
	if('annotation' in parsed):
		genes = parsed['annotation']
		# print("Number of matches in genes: "+str(compare_gene(genes,frag_select)))
		print(str(compare_gene(genes,frag_select)),end='\t')

	## test case, if PstI, compare the fragments and genes
	# if(enzyme == "PstI"):
	# 	print(compare_gene(genes))


	if(enzyme == "MspI"):
		output = open("ecoli2.fastq", "w+")
		for i in range(0,len(frag_select)):
			output.write("@Frag_"+str(i+1)+"_"+str(frag_select[i][1]+1)+"_"+str(frag_select[i][1]+150+1)+"\n")
			output.write(frag_select[i][0][:150])
			output.write("\n+\n")
			for j in range(0,150):
				output.write("F")
			output.write("\n")

			output.write("@Frag_"+str(i+1)+"_"+str(frag_select[i][2]-150+1)+"_"+str(frag_select[i][2]+1)+"\n")
			output.write(frag_select[i][0][-150:])
			output.write("\n+\n")	
			for j in range(0,150):
				output.write("F")
			output.write("\n")
		output.close()

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
	# print("Percent coverage in genome: ",end='')
	# print(count_percent_unique(genome,frag_select))
	return


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='RADSeq python script')
	parser.add_argument('-i', nargs='?', help='input genome sequence file (FASTA)')
	parser.add_argument('-db', nargs='?', default='re_db.txt',  help='resteriction enzyme dabatase file. Format per line: SbfI,CCTGCA|GG')
	parser.add_argument('-re', nargs='?', default='', help='file of list of restriction enzyme to be tested')
	parser.add_argument('-a', nargs='?', default='', help='gene annotation file for genome')
	parser.add_argument('-min', nargs='?', default=200, help='minimum fragment size (default 200)')
	parser.add_argument('-max', nargs='?', default=300, help='maximum fragment size (default 300)')
	args = parser.parse_args()

	start_time = time.time()

	minsize = args.min
	maxsize = args.max

	input_RE = ""
	parsed = {}

	## catch errors for invalid input file argument
	if args.i == None or len(args.i)==0:
		print("Sequence file not provided")
		raise SystemExit
	try:
		input_i  = open(args.i, "r+")
	except (OSError, IOError) as e:
		print("Sequence file is invalid or not found")
		raise SystemExit

	## catch errors for invalid RE DB file argument
	try:
		input_DB  = open(args.db, "r+")
		parsed['db'] = parse_enzymedb(args.db)
	except (OSError, IOError) as e:
		print("Restriction enzyme database file is invalid or not found")
		raise SystemExit
	
	## if there are annotations given by user, use it
	if args.a != None and len(args.a)>0:
		try:
			input_a  = open(args.a, "r+")
			parsed['annotation'] = parse_gff(args.a)
		except (OSError, IOError) as e:
			print("Annotation file is invalid or not found")
			raise SystemExit
	
	## if no list of preferred REs to be tested, use everything in the database, RUN HERE
	if args.re != None and len(args.re)>0:
		try:
			input_RE  = open(args.re, "r+")
			print("Name \tRE sites\tFrags filtered \tMatches in gene",end='')
			for line in input_RE:
				enz = line.strip()
				run_RE(enz, parsed, args)
		except (OSError, IOError) as e:
			print("Restriction enzyme list file is invalid or not found")
			raise SystemExit
	else:
		print("Name \tRE sites \tFrags filtered \tMatches in gene",end='')
		for key in sorted(parsed['db'].keys()):
			run_RE(key, parsed, args)

	print("\n\n--- %s seconds ---" % (time.time() - start_time))
	