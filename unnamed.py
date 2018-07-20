#!/usr/bin/python3

"""
RADSeq python script
A python script that simulates restriction enzyme digestion (single and double), size selection annd evaluates the number of fragments inside a target annotated region for given genomes in a single fasta file or generated DNA sequence with GC Frequency.
"""


from __future__ import print_function
import argparse, re, copy, operator, importlib
import sys, time, os, shutil, subprocess
import numpy as np
import remove_repeat2
from multiprocessing import Pool, Process
from multiprocessing.pool import ThreadPool

def gen(x,p):
	"""
		function that simply uses the join function on x and p (x.join(p))
	"""
	return x.join(p)

def generate_gc(gc_freq, genome_size):
	"""
		generates DNA sequence given the GC content (0<=gc<=1) and genome size.
		The generated sequence is written in a file genome.txt.
	"""
	at_freq = 1 - gc_freq
	nucleo = ['A', 'C', 'T', 'G']
	weights = [at_freq/2, gc_freq/2, at_freq/2, gc_freq/2]

	pool = ThreadPool(32)	## divide generation of bases into 4 processes
	p1 = pool.apply_async(np.random.choice, (nucleo,genome_size/4,True,weights))
	p2 = pool.apply_async(np.random.choice, (nucleo,genome_size/4,True,weights))
	p3 = pool.apply_async(np.random.choice, (nucleo,genome_size/4,True,weights))
	p4 = pool.apply_async(np.random.choice, (nucleo,genome_size-3*(genome_size/4),True,weights))
	pool.close()

	pool = ThreadPool(32)	## flatten the bases into one string
	p_1 = pool.apply_async(gen, ('',p1.get()))
	p_2 = pool.apply_async(gen, ('',p2.get()))
	p_3 = pool.apply_async(gen, ('',p3.get()))
	p_4 = pool.apply_async(gen, ('',p4.get()))
	pool.close()

	pool = ThreadPool(32)	## combine strings
	p_12 = pool.apply_async(gen, ('',[p_1.get(),p_2.get()]))
	p_34 = pool.apply_async(gen, ('',[p_3.get(),p_4.get()]))
	pool.close()

	p_1234 = ''.join([p_12.get(),p_34.get()])

	output = open("genome.txt", "w+")
	output.write(p_1234)	## write the strings into a file
	output.close()


def digest(genome, p5, p3):
	"""
		simulates the digestion process of the restriction enzymes.
		It cuts the given genome sequence with the also given p5 and p3 sites then returns a list of fragments.
	"""
	fragments = re.split(p5+p3, genome)		## split the genome into fragments

	curr_len = 0				## temporary holder for current base in genome
	new_fragments = []			## temporary list holder for digested fragments
	## Adding the hanging part
	for i in range(0,len(fragments)-1):
		temp_frag = []			## temporary list holder for current fragment
		temp_start = curr_len	## take note of curr_len, it will be the start location of the fragment
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


def shear_frag(fragments, shear_len):
	"""
		function that shears the single digest fragments. It truncates the fragments into a given shear_len.
	"""
	sheared_fragments = []
	for frag in fragments:
		temp_shear = []
		temp_frag = frag[0]
		temp_end = frag[2]
		frag_size = frag[2] - frag[1] + 1
		if (frag_size > shear_len):
			temp_frag = frag[0][:shear_len]
			temp_end = frag[2] - (frag_size - shear_len)
		temp_shear.append(temp_frag)	## fragment sequence
		temp_shear.append(frag[1])		## fragment start
		temp_shear.append(temp_end)		## fragment end
		sheared_fragments.append(temp_shear)	## append sheared fragment
	#shear_frag =[frag[:500] for frag in fragments]
	#print(sheared_fragments)
	return sheared_fragments


def dd_digest(genome_frag, p5_2, p3_2, p5, p3):
	"""
		function that simulate the double digestion.
		Returns double digested fragments
	"""
	dd_sites = 0
	dd_fragments = [];
	for i in range(0,len(genome_frag)):
		dd_frag = digest(genome_frag[i], p5_2, p3_2)
		dd_sites += len(dd_frag)
		dd_fragments.extend(dd_frag)

	## filter fragments AB+BA
	dd_filt_fragments = []
	for frag in dd_fragments:
		if (frag[0].startswith(p3_2) and frag[0].endswith(p5)) or (frag[0].startswith(p3) and frag[0].endswith(p5_2)):
			dd_filt_fragments.append(frag)

	return dd_filt_fragments

def select_size(fragments, minsize, maxsize, protocol):
	"""
		filters the fragments according to size.
		Returns a list of fragments within the given size range

	"""
	selected_fragments = []
	for frag in fragments:
		if protocol == 'ddrad':
			## fragment's start and end must be inside the gene region
			if(len(frag[0]) < maxsize and len(frag[0]) > minsize):
				selected_fragments.append(frag)
		else:
			if(len(frag[0]) > minsize):
				selected_fragments.append(frag)

	return selected_fragments


def parse_enzymedb(enzyme_db_file):
	"""
	function to parse the enzyme database file.
	Returns a dictionary of restriction enzyme and its site.
	input format: each enzyme with its restriction site in separate lines
		ex.	SbfI,CCTGCA|GG
			ApeKI,G|CWGC
	"""
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

	# if there enzymes parsed in the database file, raise an error
	if(len(list_enzymes) < 0):
		print("No restriction enzymes found in "+enzyme_db_file)
		raise SystemExit

	return list_enzymes


def parse_REinput(input_RE_file):
	"""
	function to parse the enzyme input file. Returns list of restriction enzyme's names
	format:	each enzyme in separate lines
		ex.	SbfI
			ApeKI
	"""
	input_RE  = open(input_RE_file, "r+")
	REs = []
	for line in input_RE:
		enz = line.strip()
		REs.append(enz)
	return REs

def restriction_sites(enzyme, list_enzymes):
	"""
		function that parses the restriction site. Gets the RE site given the RE's name and replaces any wildcard base.
		Returns RE sites p5 and p3 in regex format
	"""
	try:	## test if the RE is in loaded DB, raise an error if not in loaded DB
		match_enzyme = list_enzymes[enzyme]
	except:
		print("Restriction enzyme "+enzyme+" is not in database")
		raise SystemExit

	## replace wildcard bases
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


def parse_gff(annotation_file, target):
	"""
		parses the gene annotation file (gff).
		Returns list of gene locations. Per element: [0] start of gene location [1] end of gene location
	"""
	input_gff = open(annotation_file, "r+")

	gene_location = []
	for line in input_gff:
		parsed_line=line.strip().rstrip()	## strip whitespaces
		if(len(parsed_line) > 0):
			parsed_line=parsed_line.split("\t")
			if(len(parsed_line) > 2 and parsed_line[0] != '#'):
				## get only genes for now
				if(parsed_line[2] == target):			
					gene_location.append([int(parsed_line[3]),int(parsed_line[4])])	## add to gene_location the start and end
	## exit if there are no genes parsed
	if(len(gene_location) == 0):
		print("Check the annotation file. Must be in GFF format or contains none of the features wanted.")
		raise SystemExit
	return gene_location


def compare_gene(gene_location, fragments):
	"""
		counts how many fragments are within the gene region
	"""
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

def parse_input(input_name):
	"""
		function that parses the input sequence/genome file (fasta format).
		Returns a list of [0] sequence name and its corresponding [1] genome/dna sequence
	"""
	input_file  = open(input_name, "r+")

	new_genome_name = ""
	new_genome_seq = ""
	list_genome = []
	for line in input_file:
		line = line.strip().rstrip()
		if(len(line) == 0):
			continue
		if(">" in line):	## append to list after seeing ">"
			list_genome.append([new_genome_name,new_genome_seq])
			new_genome_seq =""
			new_genome_name = line
		else:
			new_genome_seq += line.strip().rstrip()	## strip whitespaces
	list_genome.append([new_genome_name,new_genome_seq])
	return list_genome[1:]


def run_RE(enzyme, parsed, args, genome):
	"""
		function that runs over the REs given a genome sequence. Performs the RADSeq process per RE
	"""
	results = open("output/"+enzyme+".out", "w+")
	results.write(enzyme+"\t")
	## if double digest
	if args.p == 'ddrad':
		enzyme1, enzyme2 = enzyme.split()
		p5,p3 = restriction_sites(enzyme1,parsed['db'])
	else:
		p5,p3 = restriction_sites(enzyme,parsed['db'])

	fragments = digest(genome, p5, p3)
	results.write(str(len(fragments))+"\t")
	## if double digest
	frag_select = []
	if args.p == 'ddrad':
		p5_2, p3_2 = restriction_sites(enzyme2,parsed['db'])
		dig_frag = [item[0] for item in fragments]
		fragments = dd_digest(dig_frag,p5_2,p3_2,p5,p3)
		frag_select = select_size(fragments,int(args.min), int(args.max),args.p)

	## if single digest, shear fragments
	else:
		frag_select = select_size(fragments,int(args.min), int(args.max),args.p)
		frag_select = shear_frag(frag_select,int(args.max))

	results.write(str(len(frag_select))+"\t")
	# ## select the fragments based on size
	# frag_select = list(filter(lambda frag: (frag[2]-frag[1]) < maxsize and (frag[2]-frag[1]) > minsize, shear_frag))
	#frag_select = select_size(shear_frag, args.min, args.max)

	

	output = open("reads/"+enzyme+"_read1.fastq", "w+")
	output2 = open("reads/"+enzyme+"_read2.fastq", "w+")
	for i in range(0,len(frag_select)):
		output.write("@Frag_"+str(i+1)+"_"+str(frag_select[i][1]+1)+"_"+str(frag_select[i][2]+1)+"\n")
		output.write(frag_select[i][0][:100])
		output.write("\n+\n")
		for j in range(0,100):
			output.write("A")
		output.write("\n")

		output2.write("@Frag_"+str(i+1)+"_"+str(frag_select[i][1]+1)+"_"+str(frag_select[i][2]+1)+"\n")
		output2.write(frag_select[i][0][-100:])
		output2.write("\n+\n")
		for j in range(0,100):
			output2.write("A")
		output2.write("\n")
	output.close()
	output2.close()

	shellscript = subprocess.Popen(["./bwa_aln.sh %s %s" % (args.i,enzyme)], shell=True, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, close_fds=True)
	shellscript.wait()

	unique_repeats,uniq_count,rept_count = remove_repeat2.remove_XAs(enzyme)

	results.write(str(uniq_count)+"\t"+str(rept_count)+"\t"+str(unique_repeats)+"\t")
	if('annotation' in parsed):
		genes = parsed['annotation']
		# print("Number of matches in genes: "+str(compare_gene(genes,frag_select)))
		hit_genes = compare_gene(genes,frag_select)
		print(enzyme+'\t'+str(len(fragments))+'\t'+str(len(frag_select))+"\t"+str(hit_genes))
		results.write(str(hit_genes)+"\t")
		results.write(str())
	else:
		print(enzyme+'\t'+str(len(fragments))+'\t'+str(len(frag_select))+"\t"+str(uniq_count)+"\t"+str(rept_count)+"\t"+str(unique_repeats))

	results.write("\n")
	results.close()
	return


def run_genome(REs, parsed, args,list_genomes):
	"""
		function that calls the run_RE function over multiple genomes/sequences
	"""


	for i in range(0,len(list_genomes)):
		genome = list_genomes[i][1]
		print("FASTA: "+list_genomes[i][0][1:])
		if('annotation' in parsed):
			print("Name \tRE sites\tFrags filtered \tMatches in gene")
		else:
			print("Name \tRE sites\tFrags filtered")
		pool = Pool(32)
		for enz in REs:
			# run_RE(enz,parsed,args,genome)
			p = Process(target=run_RE,args=(enz,parsed,args,genome))
			p.start()
			p.join()

if __name__ == '__main__':

	## argument parser
	parser = argparse.ArgumentParser(description='RADSeq python script')
	parser.add_argument('-i', nargs='?', help='input genome sequence file (FASTA)')
	parser.add_argument('-db', nargs='?', default='re_db.txt',  help='resteriction enzyme dabatase file. Format per line: SbfI,CCTGCA|GG')
	parser.add_argument('-re', nargs='?', default='', help='file of list of restriction enzyme to be tested')
	parser.add_argument('-a', nargs='?', default='', help='gene annotation file for genome')
	parser.add_argument('-at', nargs='?', default='gene', help='what to look for in gene annotation file (ex. gene region, exon, intron, etc)')
	parser.add_argument('-min', nargs='?', default=200, help='minimum fragment size (default 200)')
	parser.add_argument('-max', nargs='?', default=300, help='maximum fragment size (default 300)')
	parser.add_argument('-p', nargs='?', default='orig', help='radseq protocol: use ddrad for double digestion')
	parser.add_argument('-gc', nargs='?', help='input gc frequency. Value must be between 0 and 1')
	parser.add_argument('-dna', nargs='?', help='input dna estimated length')

	args = parser.parse_args()

	start_time = time.time()

	input_RE = ""
	parsed = {}
	genome = ""

	importlib.import_module("remove_repeat2")
	if(os.path.exists("reads") == True):
		shutil.rmtree("reads")
	os.makedirs("reads")

	if(os.path.exists("output") == True):
		shutil.rmtree("output")
	os.makedirs("output")


	## catch errors for invalid input file argument
	if (args.i == None and args.gc != None and args.dna != None):
		print("Choose only one input source")
		raise SystemExit
	if (args.i == None or len(args.i)==0) and (args.gc == None and args.dna == None):
		print("Sequence file / Input not provided")
		raise SystemExit

	## if unknown genome and given gc frequency, generate genome sequence
	if (args.gc != None and args.dna != None) and args.i == None:
		gc_freq = float(args.gc)
		if(gc_freq <= 1 and gc_freq >= 0):
			genome_length = int(args.dna)
			generate_gc(gc_freq,genome_length)
			genome = parse_input('genome.txt')
		else:
			print("GC frequency must be between 0 and 1")
			raise SystemExit

	## no given gc frequency or has input genome, open genome file
	else:
		try:
			input_i  = open(args.i, "r+")
			input_i.close()
			shellscript = subprocess.Popen(["./bwa_index.sh %s" % args.i], shell=True, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, close_fds=True)
			genome = parse_input(args.i)
			
			
			# subprocess.call("./bwa_index.sh %s" % args.i,shell=)
		except (OSError, IOError) as e:
			print("Sequence file "+args.i+" is invalid or not found")
			raise SystemExit


	## try to open the RE DB file, catch errors found
	try:
		input_DB  = open(args.db, "r+")
		parsed['db'] = parse_enzymedb(args.db)
	except (OSError, IOError) as e:
		print("Restriction enzyme database file is invalid or not found")
		raise SystemExit

	## catch errors for invalid protocol
	if (args.p != 'orig' and args.p != 'ddrad'):
		print("Invalid RADSeq protocol. Use 'orig' for Original RADSeq, 'ddrad' for ddRADSeq")
		raise SystemExit

	## require RE file if protocol is ddrad
	if (args.p == 'ddrad'):
		if(args.re == None or len(args.re) < 1):
			print("DDRad Protocol requires an -re argument")
			raise SystemExit

	## if there are annotations given by user, use it (parse it)
	if (args.a != None and len(args.a)>0):
		try:
			input_a  = open(args.a, "r+")
			parsed['annotation'] = parse_gff(args.a, args.at)
		except (OSError, IOError) as e:
			print("Annotation file is invalid or not found")
			raise SystemExit


	##### running the core of the script happens here #####
	shellscript.wait()
	# if protocol is original and RE file is given
	# parse the RE file then call run_genome
	if (args.re != None and len(args.re) > 0 and args.p == 'orig'):
		try:
			REs  = parse_REinput(args.re)
			run_genome(REs, parsed, args,genome)
		except (OSError, IOError) as e:
			print("Restriction enzyme list file is invalid or not found")
			raise SystemExit

	## if protocol is ddrad and RE file is given
	## parse the RE file then call run_genome
	elif (args.re != None and len(args.re) > 0 and args.p == 'ddrad'):
		try:
			REs  = parse_REinput(args.re)
			run_genome(REs, parsed, args,genome)
		except (OSError, IOError) as e:
			print("Restriction enzyme list file is invalid or not found")
			raise SystemExit

	## if no RE file given, use everything in the database and protocol is original by default
	else:
		run_genome(sorted(parsed['db'].keys()), parsed, args,genome)

	print("\n\n--- %s seconds ---" % (time.time() - start_time))
