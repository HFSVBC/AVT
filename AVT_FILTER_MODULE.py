#!/usr/local/bin/python
# coding: utf-8
"""
Universidade de Lisboa
Faculdade de Ciências
Departamento de Informática
LaSIGE

Allele Validation Tool (AVT) - FILTER MODULE

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

--------------------------------------------------------------------------------------------
EXECUTION METHOD:

Exemple:
./AVT_FILTER_MODULE.py <AlleleCallresults.txt> <ContigLocation.txt> <file.fasta> <out.fasta>
"""

__author__ = "Hugo Filipe Curado, Margarida Cândido"
__copyright__ = "2016 under a "
__version__ = "1.0"
__maintainer__ = "Hugo Filipe Curado"
__email__ = "hugofsvbc@alunos.fc.ul.pt"

import sys
import re
from copy import deepcopy

#funcao para abrir ficheiros e separar as duas linhas
def FileOpen (nome):
	"""
	Opens txt files with only two lines

	Requires: the name of the txt file as a string, and that the file only has two lines.
	Ensures: two lists with each line of the file separated by ' '.
	"""
	with open(nome, 'r') as fp:
		reader = (fp.read()).split('\n')
		#separa as duas linhas
		line1 = (reader[0]).split('	')
		line2 = (reader[1]).split('	')
	#retorna as duas linhas
	return line1, line2

def reference_grabber(AlleleCallresults):
	"""
	Gets the references from the given file that have an 'INF-' on their colun.

	Requires: The name of valid AlleleCallresults.txt as string.
	Ensures: A list with all the references corresponding to the ones that have
	an 'INF-' on their collun.
	"""

	#the output of FileOpen, runned with the first file
	line1, line2 = FileOpen(AlleleCallresults)
	i = 0
	index = []

	for value in line2:
		if re.search("INF-", value):
			index=index+[i]
		i+=1

	references = []

	for ind in index:
		references = references + [line1[ind]]

	return references

def NODE_grabber(ContigLocation, references_alleles):
	"""
	From a list of references (references_alleles) finds the NODEs in the same colun.

	Requires: The name of the txt ContigLocation file (exemple/ContigLocation.txt), as a
	string, and a list with references of the alleles references_alleles.
	Ensures: A dictionary with the references as key and a list with the NODE, biggining
	position and end position of the bit of sequence to get, as a value.
	"""
	line1, line2 = FileOpen(ContigLocation)
	i = 0
	index = []
	for value in line1:
		for ref in references_alleles:
			if value == ref:
				index = index + [[ref, i]]
		i+=1

	NODEs = {}
	for li in index:
		NODEs [li[0]] = line2[(li[1])].split('&')
	for key in NODEs:
		NODEs[key] = [NODEs[key][0], (NODEs[key][1]).split('-')]
	for key in NODEs:
		NODEs[key] = [NODEs[key][0], int(NODEs[key][1][0])-200, int(NODEs[key][1][1])+201] 

	return NODEs

def filt_seq_grabber(file_seq, references_node):
	"""
	From the NODEs in a dictionary extracts from the .fasta file the bit of sequence to study.

	Requires: The name of a txt file_seq file (exemple/10_S10_L001.fasta), as a str, and a
	dictionatry with the references of the alleles as key and as value a list with the NODEs to
	get as well as the biggining and end position of the bit of sequence to get.
	Ensures: A dictionatry with the reference as key and an str as value with the bit of sequence
	to study.
	"""
	Seq = {}
	key = []
	value = []
	with open(file_seq, 'r') as fp:
		reader = (fp.read()).split('\n')
		for line in reader:
			if re.search(">", line):
				key = line
				Seq[key] = value
				key = ''
				value = ''
			else:
				value = value + line.replace("\n", "")
	out={}
	for key in Seq:
		for ref in references_node:
			if key == '>'+references_node[ref][0]:
				out[ref] = Seq[key][(references_node[ref][1]):(references_node[ref][2])]
	return out


def FileWriter(name, output):
	"""
	From a dictionary writes a .fasta file with a given name (exemple/10_S10_L001_filt.fasta).

	Requires: The name of the .fasta file to write, and a dictionary with the references of the
	alleles as a key and the bits of sequence to analyse as a value
	Ensures: A .fasta file with the given name conatining the sequence to analyse and the reference.
	"""
	with open(name, 'w') as fp:
		i=0
		out = ''
		for ref in output:
			out=''
			i=0
			for value in output[ref]:
				if i==60:
					out = out + '\n'
					i=0
				out = deepcopy(out + value)
				i+=1
			output[ref]=out
			fp.write('>'+ref+'\n'+str(output[ref])+'\n')
	fp.close()
	

def filt(file1, file2, file3, nameOut):
	"""
	Separates the bits of sequence to analyse from the whole sequence

	Requires: the name of the three files as strings, the first two as txt and the third as fasta.
	A string for the name of the file to output.
	Ensures: a file with the given name with the bits of sequence to analyse (exemple: exemple/
	10_S10_L001_filt.fasta.
	"""
#29041996
	references = reference_grabber(file1)
	NODES = NODE_grabber(file2, references)
	filt_seq = filt_seq_grabber(file3, NODES)
	FileWriter(nameOut, filt_seq)