"""
RApyDS
Restriction Site Associated DNA Python-Digested Simulation

to_json.py
"""
import io, json

def convert_json(genome_name):
	whole_json = {}
	cut_json={}

	with open("output/csv_"+genome_name+".csv") as f:
		for content in f:
			line = content.split("\t")
			temp_arr = []
			## create lists of fragment's length and cut sites
			whole_json[line[0]] = [int(x) for x in line[1].split(",")]
			cut_json[line[0]] = [int(x) for x in line[2].split(",")]
	
	## dump fragment length list to json file
	with open("output/"+genome_name+'.json', 'w') as fp:
		json.dump(whole_json, fp)

	## dump cut site list to json file
	with open("output/cut"+genome_name+'.json', 'w') as fp:
		json.dump(cut_json, fp)

if __name__ == '__main__':
	pass