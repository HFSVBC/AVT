#!/usr/local/bin/python
# coding: utf-8
"""
Universidade de Lisboa
Faculdade de Ciências
Departamento de Informática
LaSIGE

Allele Validation Tool (AVT) - PREPARATION MODULE
CHANGE THIS
Script capable of separating the alleles to validate from DNA sequence in fasta
file (exemple: exemple/10_S10_L001.fasta).
The script needs three inputs. A txt file with all the inf and references (exemple:
exemple/AlleleCallresults.txt). A second txt file with the nodes and it's refrences
(exemple: exemple/ContigLocation.txt). And a third file with the sequence (exemple:
exemple/10_S10_L001.fasta). The output will be a fasta file with the bits of sequence
to analyse with the appropriate threshold.

--------------------------------------------------------------------------------------------
SYSTEM REQUIREMENTS:

-UNIX based OS or Windows
-BWA installed in the /usr/local/bin directory
-SAMTOOLS installed in the /usr/local/bin directory
-BCFTOOLS installed in the /usr/local/bin directory

"""

__author__ = "Hugo Filipe Curado, Margarida Cândido"
__copyright__ = ""
__version__ = "1.0"
__maintainer__ = "Hugo Filipe Curado"
__email__ = "hugofsvbc@alunos.fc.ul.pt"

import os
from subprocess import call
from AVT_FILTER_MODULE import filt

path = ''
path_out = ''

def avtFiltration(AlleleCallresults, ContigLocation, RefFasta):
	"""
	This function calls the de novo assembler result filtration module

	Requires: the name of three files as string (AlleleCallresults, ContigLocation,
	RefFasta).
	"""
	global path, path_out
	path = os.path.join('genData', 'ref_filt.fasta')
	print "--*--REFERENCE SEQUENCE FILTRATION STARTED--*--"
	filt(AlleleCallresults, ContigLocation, RefFasta, path)
	print "--*--REFERENCE SEQUENCE FILTRATION Ended--*--\n"

def avtBWAmap(Reads1, Reads2):
	"""
	This function calls the BWA Mapper to map the output of the avtFiltration
	against the reads, if paired ended, or the read, if single ended.

	Requires: the name of the file(s) of the read(s) as string
	Ensures: a typed sam file with the result of the map.
	"""
	global path, path_out
	print "--*--BWA MAPPING STARTED--*--"
	call(["bwa", "index", path])
	path_out = os.path.join('genData', 'aln-pe.sam')
	with open(path_out, 'w') as f:
		print "\n"
		#--------------VERIFIES IF THE -R2 OPTION WAS SET
		if Reads2 != '':
			call(["bwa", "mem", path, Reads1, Reads2], stdout=f)
		else:
			call(["bwa", "mem", path, Reads1], stdout=f)

	print "--*--BWA MAPPING ENDED--*--\n"

def avtBAMconversiom():
	"""
	This function converts the sam file into a binary file, BAM

	Ensures: a bam file
	"""
	global path, path_out
	print "--*--BAM CONVERSION STARTED--*--"
	path = path_out
	path_out = os.path.join('genData', 'aln-pe.bam')
	with open(path_out, 'wb') as f:
		call(["samtools", "view", "-Sb", path], stdout=f)
	path = path_out
	path_out = os.path.join('genData', 'aln-pe_SORTED')
	call(["samtools", "sort", path, path_out])
	path = os.path.join('genData', 'aln-pe_SORTED.bam')
	call(["samtools", "index", path])
	print "--*--BAM CONVERSION ENDED--*--\n"

def avtVCFconversion():
	"""
	Converts the bam file into a VCF, variant call format, file

	Ensures: A VCF file
	"""
	global path, path_out
	print "--*--VCF CONVERSION STARTED--*--"
	path_out = os.path.join('genData', 'aln-pe_SORTED.bcf')
	with open(path_out, 'w') as f:
		call(["samtools", "mpileup", "-g", "-f", os.path.join('genData', 'ref_filt.fasta'), path], stdout=f)
	path = path_out
	path_out = os.path.join('genData', 'aln-pe_SORTED.vcf')
	with open(path_out, 'w') as f:
		call(["bcftools", "call", "-c", "-v", path], stdout=f)
	print "--*--VCF CONVERSION ENDED--*--\n"
