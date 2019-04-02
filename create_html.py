"""
RApyDS
Restriction Site Associated DNA Python-Digested Simulation

create_html.py
"""

from __future__ import print_function
import shutil, os, errno

footer = "\n\t</body>\n</html>"
MAX = 5

def create_overview():
	"""
		function that creates the overview.html file
	"""
	file_html = open("index.html", "w+")
	
	with open("templates/head.txt", "r") as fp:
		for content in fp:
			file_html.write(content)

	genome_names = []
	with open("output/genome_names.txt", "r") as f:
		for content in f:
			genome_names.append(content.strip())

	file_html.write('<div class="three wide column"><div class="ui vertical fluid tabular menu">')
	side_bar = '\t\t\t\t<a class="item%s" data-tab="%s">%s</a>\n'
	active = False
	for name in genome_names:
		file_name = name.split(" ")[0]
		if(active == True):
			file_html.write(side_bar % ("",file_name, file_name))
		else:
			file_html.write(side_bar % (" active",file_name, file_name))
			active = True
	active = False

	file_html.write("\t\t\t</div>\n\t\t</div>\n")
	file_html.write('\n\t\t\t<div class="thirteen wide stretched column">\n\t\t\t\t<div class="ui segment">\n')

	header_name = '\t\t\t\t\t<div class="ui tab%s" data-tab="%s"><h2>%s</h2>\n'
	table_header = '<table id="" class="ui celled table enzymes_table" style="width:100%">\n\t<thead>\n<tr>\n\t\t<th>Enzyme</th>\n\t\t<th>Fragments after Digestion</th>\n\t\t<th>Fragments after Shearing and Size Select</th>\n\t\t<th>Percent Coverage</th>\n\t\t<th>\n<div data-tooltip="Reads that no repeats are found after running BWA">Unique Reads <i class="icon small question circle"></i></div></th>\n\t\t<th><div data-tooltip="Reads with repeats that are found after running BWA">Repeat Reads <i class="icon small question circle"></i></div></th>\n\t\t<th><div data-tooltip="Regions found by BWA that are repeats of reads">Repeat regions hit by Reads <i class="icon small question circle"></i></div></th>\n\t\t<th>Fragments within Annotated Region</th>\n\t\t<th>Annotated Regions hit by Fragments</th>\n\t\t<th>Percent of Annotated Covered</th></tr>\n\t</thead>\n\t<tbody>\n'
	table_def = '<td>%s</td>'
	table_def_per = '<td>%s%%</td>'

	for name in genome_names:
		file_name = name.split(" ")[0]
		if(active == True):
			file_html.write(header_name % ("", file_name, name))
		else:
			file_html.write(header_name % (" active", file_name, name))
			active = True

		file_html.write(table_header)


		with open("output/"+file_name+".txt", "r") as results_file:
			for result_line in results_file:
				results = result_line.strip().split("\t")
				file_html.write("<tr>")
				for i in range(len(results)):
					# if(i == 4 or i == 9):
					file_html.write(table_def % results[i])
				file_html.write("</tr>\n")

		file_html.write("</tbody></table>")
		file_html.write("\n\t\t\t\t\t</div>\n")
	

	with open("templates/foot.txt", "r") as fp:
		for content in fp:
			file_html.write(content)

	with open("templates/overview_js.txt", "r") as fp:
		for content in fp:
			file_html.write(content)

	file_html.write(footer)

	file_html.close()


def create_gel_html():
	"""
		function that creates the gel.html file
	"""
	file_html = open("gel.html", "w+")
	
	with open("templates/head.txt", "r") as fp:
		for content in fp:
			file_html.write(content)


	file_html.write('\n\t\t\t\t<div class="eight wide column">\t\t\t\t<form class="ui form">\n\t\t\t\t\t<div class="field">\n\t\t\t\t\t<label>Genome</label>\n\t\t\t\t\t<select id="genome">\n\t\t\t\t\t<option value="" selected=""></option>\n')

	with open("output/genome_names.txt", "r") as fp:
		genome_select = '\t\t\t\t\t<option value="%s">%s</option>\n'
		for content in fp:
			line = content.strip()
			file_html.write(genome_select % (line.split(' ')[0], line))
	file_html.write('\t\t\t\t\t</select>\n\t\t\t\t\t</div>\n\n')
	file_html.write('\t\t\t\t<div class="two fields">\n\t\t\t\t\t\t<div class="field"><label>Start</label><input type="text" class="three_fields" id="gel_start" value="0" size="3"/></div>\n\t\t\t\t\t\t<div class="field"><label>End (set both to 0 for no limit)</label><input type="text" class="three_fields"  id="gel_end" value="0" size="3"/></div></div>\n\t\t\t\t\t<div class="field"><label>Markers</label><input type="text" class="three_fields" id="gel_step" value="10,50,100,500,1000,2000,3000,4000"size="3"/></div>')
	file_html.write('\n\t\t\t\t\t<div class="five fields">')
	enzymes = []
	with open("output/RE.txt", "r") as fp:
		for content in fp:
			line = content.strip()
			enzymes.append(line)

	for i in range(0,MAX):
		file_html.write('\n\n\t\t\t\t\t<div class="field">\n\t\t\t\t\t\t<select id="enzyme_%s">\n\t\t\t\t\t\t\t<option value="" selected=""></option>' % str(i+1))
		for enz in enzymes:
			enz_option = '\n\t\t\t\t\t\t\t<option value="%s">%s</option>'
			file_html.write(enz_option % (enz, enz))
		file_html.write('\n\t\t\t\t\t\t</select>\n\t\t\t\t\t</div>')


	file_html.write('\n\t\t\t\t</div>')
	file_html.write('\n\t\t\t\t<input type="button" value="Simulate" onclick="read_data()" /><br/>\n\t\t\t\t</form>\n\t\t\t\t</div>\n\t\t\t\t<div class="eight wide column">\n\t\t\t\t\t<div class="gel_body" style="width: 600px; height: 650px; background-color: #F0F0F0;"></div></div>')
	with open("templates/foot.txt", "r") as fp:
		for content in fp:
			file_html.write(content)

	with open("templates/gel_js.txt", "r") as fp:
		for content in fp:
			file_html.write(content)

	file_html.write(footer)
	file_html.close()

def create_cutsite():
	"""
		function that creates the cutsite.html file
	"""
	shutil.copyfile('templates/cutsite.html', 'cutsite.html')


def create_report(output_name):
	"""
		function that creates the report files and zips them
	"""
	## call the functions to create the html files
	create_overview()
	create_gel_html()
	create_cutsite()

	## create output directory
	if(os.path.exists(output_name) == True):
		shutil.rmtree(output_name)
	os.makedirs(output_name)
	try:
	    os.remove(output_name+'.zip')
	except OSError:
	    pass

	## copy html files and directories src and output
	shutil.move('index.html', output_name)
	shutil.move('gel.html', output_name)
	shutil.move('cutsite.html', output_name)

	shutil.copytree('output', output_name+'/output/')
	shutil.copytree('src', output_name+'/src')
	# try:
	# 	shutil.copytree('output', output_name+'/output/')
	# except:
	# 	print
	# try:
	# 	shutil.copytree('src', output_name+'/src')
	# except:
	# 	pass

	## make the zip file then do cleanup
	
	shutil.make_archive(output_name, 'zip', output_name)
	if(os.path.exists(output_name) == True):
		shutil.rmtree(output_name)


if __name__ == '__main__':
	create_report("report")