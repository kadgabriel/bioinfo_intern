
## Restriction Site Associated DNA Python-Digested Simulation  (RApyDS)

> still under development

**RApyDS** is a python script that performs in silico digestion as an aid for choosing a restriction enzyme for RADseq. 

**RApyDS** can perform in silico digestion of a given genome or known gc content. It can simulate single-enzyme or double-enzyme digestion on a single file or multiple FASTA files. After digestion, the program can match fragments to a given gene annotation file and run alignment analysis to check for unique and repeat regions.

**RApyDS** provides a detailed report and visualization of the simulated digestion. 


### Documentation
- For the [RApyDS program overview / documentation](docs/rapyds.pdf)
- For the [full RApyDS program technical manual](docs/rapyds_manual.pdf)

### Requirements

- Python 2.7 or greater
- numpy ( `pip install numpy` )
- BWA 0.7.12 (http://bio-bwa.sourceforge.net/)
- Firefox (for viewing the html files)
- Linux OS

### Usage

#### Arguments
-  **-h, --help**  show this help message and exit
-  **-i [I]**      input genome sequence file (FASTA)
-  **-db [DB]**    resteriction enzyme dabatase file. Format per line: SbfI,CCTGCA|GG
-  **-re [RE]**    file of list of restriction enzyme to be tested
-  **-a [A]**      annotation file for genome (GFF)
-  **-at [AT]**    what to look for in gene annotation file (ex. gene region, exon, intron, etc) (default: gene)
-  **-min [MIN]**  minimum fragment size (default: 200)
-  **-max [MAX]**  maximum fragment size (default: 300)
-  **-bp [BP]**    base pair read length for FASTQ generation (default: 100)
-  **-p [P]**      radseq protocol: use ddrad for double digestion (default: orig)
-  **-gc [GC]**    input gc frequency. Value must be between 0 and 1
-  **-dna [DNA]**  input dna estimated length
-  **-o [O]**      output file name (default: report)
-  **-t [T]**      number of processes (default 4)



#### Common Usage
- To use RApyDS with a given genome file (and/or corresponding annotation file)

``python rapyds.py -i <input file.fasta> [-a <annotation file.gff>] [other args]``

- To use RApyDS by generating a genome file given GC content/frequency

``python rapyds.py -gc <gc frequency from 0 to 1> -dna <length> [other args]``


#### Sample Run
Given an E.Coli FASTA ``ecoli.fasta`` with annotation file ``ecoli_annotation.gff``, the RADSeq protocol is DDRad

``python rapyds.py -i ecoli.fasta -a ecoli_annotation.gff -p ddrad``


#### Output

The output files of the program is a zip file named after the ``-o`` argument or by default, the ``report.zip``.
Inside the archive are 3 html files containing:
- Overview (``index.html``) - summary of the fragments' data
- Electrophoresis (``gel.html``) - electrophoresis simulation comparing a genome and choice of up to 5 restriction enzyme
- Cut Site Distribution (``cutsite.html``) - zoomable images of the genome marked with cut sites by choice of enzymes

> It is advisable open the html files using Firefox. There is an issue with Google Chrome when opening local files. Only the overview file will work fine for any browser.

> Slowdown may be experienced when loading visualisations for enzymes with large number of cut sites.

See sample output [here](docs/examples/)

### Authors

**RADSeq Team** - Arielle Gabriel. Mark Mendoza. Danielle Pamulaklakin

Project for the 2018 Internship Program in Bioinformatics