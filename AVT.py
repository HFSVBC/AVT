#!/usr/local/bin/python
# coding: utf-8
"""
Universidade de Lisboa
Faculdade de Ciências
Departamento de Informática
LaSIGE

This file is part of AVT.

AVT is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

AVT is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with AVT.  If not, see <http://www.gnu.org/licenses/>.


Allele Validation Tool (AVT) - MAIN MODULE

Script capable of separating the alleles to validate from DNA sequence in fasta
file (example: example/10_S10_L001.fasta).
The script needs four / five inputs. A txt file with all the inf and references (example:
example/AlleleCallresults.txt). A second txt file with the nodes and it's refrences
(example: example/ContigLocation.txt). A third file with the sequence (example:
example/10_S10_L001.fasta). And the read(s) (example: example/10_S10_L001_R1_001/
and example/10_S10_L001_R2_001/). The output will be a txt file with
the alleles that have not passed the validation."

--------------------------------------------------------------------------------------------
SYSTEM REQUIREMENTS:

-UNIX based OS or Windows
-BWA installed in the /usr/local/bin directory
-SAMTOOLS installed in the /usr/local/bin directory
-BCFTOOLS installed in the /usr/local/bin directory

--------------------------------------------------------------------------------------------
EXECUTION METHOD:
./AVT <AlleleCallresults.txt> <ContigLocation.txt> <file.fasta> <reads.fq> <options>

or

./AVT <AlleleCallresults.txt> <ContigLocation.txt> <file.fasta> <read1.fq> -R2 <read2.fq> <options>
example:

"""

__author__ = "Hugo Filipe Curado, Margarida Cândido"
__copyright__ = """
	This file is part of AVT.

    AVT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AVT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AVT.  If not, see <http://www.gnu.org/licenses/>."""
__version__ = "1.0"
__maintainer__ = "Hugo Filipe Curado"
__email__ = "hugofsvbc@alunos.fc.ul.pt"

import os
import AVT_PREPARATION_MODULE as AVTpm

#----------------AVT's Core functions----------------
#----------------AVT's Auxiliar Core functions----------------
def avtVCFdataprep(reader):
	"""
	This function jumps all of the lines that start with a # on a VCF file
	and stores the probality of error of an Alt in the VCF file in a dictionary

	Requires: a reader of a well strutured VCF file
	Ensures: a dictionary with the probality of error of an Alt in the VCF file
	in a dictionary. The key is a tuple with the id of the allele and the position
	of the SNP.
	"""
	i=0
	while ((reader[i]).startswith('#')):
		i+=1

	reader = reader[i:]
	
	out = {}
	for x in reader:
		x = x.split('\t')
		if x != ['']:
			out[(x[0],x[1])]=x[2:6]
			out[(x[0],x[1])][3]=10**(-float(out[(x[0],x[1])][3])/10)

	return out

def avtVCFfilterBadAlleles(vcfTable, threshold):
	"""
	This function filters all of the probalities of error greather than the threshold.

	Requires: a dictionary with the probality of error of an Alt in the VCF file.
	Ensures: a list with all of the alelles that have a probality of error less than
	the threshold.
	"""
	return filter(lambda x: vcfTable[x][3] < threshold, vcfTable.keys())

def avtPrintBadAlleles(badAlelles, lines, output):
	"""
	This function writes to txt file the alleles that the tool thinks were
	bad calls.

	Requires: a list containing all of the alelles to write to file
	Ensures: a txt file similar to VCF format with the all the bad calls
	"""
	with open(output, 'w') as fp:
		writer = '#CHROM\tPOS\tID\tREF\tALT\tERROR_PERCENTAGE(ALT_is_Wrong)\n'
		for line in badAlelles:
			writer += line[0]+'\t'+line[1]+'\t'+("\t".join(str(x) for x in lines[line][0:3]))+'\t'+str(lines[line][3]*100)+'\n'
		fp.write(writer)

def avtVCFanalyses(threshold, output):
	"""
	This function reads the vcf file and separates the data based on QUAL values and
	a threshold

	Requires: a threshold, float, between 0 and 1
	"""
	with open(os.path.join('genData', 'aln-pe_SORTED.vcf'), 'r') as fp:
		reader = (fp.read()).split('\n')
		if reader[-1] == '':
			reader = reader[:-1]

	data = avtVCFdataprep(reader)

	badAlelles = avtVCFfilterBadAlleles(data, threshold)

	avtPrintBadAlleles(badAlelles, data, output)

#----------------AVT's Main Core function----------------
def avtMAIN(AlleleCallresults, ContigLocation, RefFasta, Reads1, Reads2, threshold, output, delete):
	"""
	The main module calls of the functions needed to get a txt file with all of the bad allignments.

	Requires: The name of the files containing the AlleleCallresults, ContigLocation, RefFasta, Reads1,
	Reads2 all in string, the threshold in float, the output file name in string and the delete as boolean.
	Ensures: A txt file, similar to the VCF format, with the name contained in output with all the bad aligments.
	"""
	#--------------CREATES FOLDER FOR GENERATED DATA
	if not os.path.exists('genData'):
		os.mkdir('genData')
	else:
		for File in os.listdir('genData'):
			os.remove(os.path.join('genData', File))
	#--------------DATA PREP
	AVTpm.avtFiltration(AlleleCallresults, ContigLocation, RefFasta)
	AVTpm.avtBWAmap(Reads1, Reads2)
	AVTpm.avtBAMconversiom()
	AVTpm.avtVCFconversion()
	print "--*--VCF ANALYSES STARTED--*--"
	avtVCFanalyses(threshold, output)
	print 'Data writen to', output
	print "--*--VCF ANALYSES ENDED--*--"
	if delete == 'True':
		for File in os.listdir('genData'):
			os.remove(os.path.join('genData', File))
		os.rmdir("genData")

#----------------Var retrival and options check----------------
if __name__ == '__main__':
    import argparse, os
    parser = argparse.ArgumentParser(description = "Allele Validation Tool (AVT) - \
		Script capable of separating the alleles to validate from DNA sequence in fasta\
		file (example: example/10_S10_L001.fasta).\
		The script needs four / five inputs. A txt file with all the inf and references (example:\
		example/AlleleCallresults.txt). A second txt file with the nodes and it's refrences\
		(example: example/ContigLocation.txt). A third file with the sequence (example:\
		example/10_S10_L001.fasta). And the read(s) (example: example/10_S10_L001_R1_001/\
		and example/10_S10_L001_R2_001/). The output will be a txt file with\
		the alleles that have not passed the validation.")
    #INPUT FILES
    #Allele_Call_results
    parser.add_argument('Allele_Call_results', metavar = 'AlleleCallresults.txt',
                        help = 'the file containing the allele call results')
    #Contig_Location
    parser.add_argument('Contig_Location', metavar = 'ContigLocation.txt',
                        help = 'the file containing the contig locatio')
    #In_Sequence
    parser.add_argument('In_Sequence', metavar = 'file.fasta',
                        help = 'the file containing the sequence generated by the de novo aligner')
    #First_Read
    parser.add_argument('First_Read', metavar = 'read1.fq',
                        help = 'the file containing the read number one')
    #OPTIONS
    #Second_Read
    parser.add_argument('-R2, --Second_Read', metavar = 'R2', dest = 'Second_Read',
                        default = "", type = str,
                        help = '<read2.fq> the file containing the read number two')
    #threshold
    parser.add_argument('-t, --threshold', metavar = 't', dest = 'threshold',
                        default = 0.1, type = float,
                        help = 'the threshold to consider the alleles as wrong calls')
    #output
    parser.add_argument('-o, --output', metavar = 'o', dest = 'output',
                        default = "output.txt", type = str,
                        help = 'the name of the output file')
    #delete_genData
    parser.add_argument('-d, --delete', metavar = 'd', dest = 'delete',
                        default = 'False', type = str,
                        help = '<True> deletes the genData folder after analyses end')

    arguments = parser.parse_args()
    avtMAIN(arguments.Allele_Call_results, arguments.Contig_Location, arguments.In_Sequence, arguments.First_Read,\
    	arguments.Second_Read, arguments.threshold, arguments.output, arguments.delete)
