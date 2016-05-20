# AVT - Allele Validation Tool (AVT) for whole genome Multilocus Sequence Typing (wgMLST) schemas
AVT (Allele Validation Tool), is a tool capable of validate novel alleles found in de novo assemblers by allele calling algorithms.

AVT's website: http://hfsvbc.github.io/AVT/

###Install
Simply downaload a [zip of the repository](https://github.com/HFSVBC/AVT/archive/master.zip) and install all of AVT's <a name="SR"></a>dependencies.

###How to Use
For single ended reads
```
python AVT.py <AlleleCallresults.txt> <ContigLocation.txt> <file.fasta> <reads.fq>
```
For paired ended reads
```
python AVT.py <AlleleCallresults.txt> <ContigLocation.txt> <file.fasta> <read1.fq> -R2 <read2.fq>
```
Available opthions:
```
-R2, --Second_Read R2   <read2.fq> the file containing the read number two
 -t, --threshold t      the threshold to consider the alleles as wrong calls
 -o, --output o         the name of the output file
 -d, --delete d         <True> deletes the genData folder after analyses end
```
For the help menu
```
python AVT.py -h
```
###System Requirements
- UNIX based system
- [BWA](http://bio-bwa.sourceforge.net) installed in the /usr/local/bin directory
- [SAMTOOLS](http://samtools.sourceforge.net) installed in the /usr/local/bin directory
- [BCFTOOLS](http://www.htslib.org/download/) installed in the /usr/local/bin directory

###Authors
Hugo Curado(1), Margarida Cândido(1), Mickael Silva(2), João André Carriço(2) e Francisco M. Couto(1)

(1) LaSIGE, Departamento de Informática, Faculdade de Ciências, Universidade de Lisboa, Lisboa, Portugal 

(2) Instituto de Microbiologia, Instituto de Medicina Molecular, Faculdade de Medicina, Universidade de Lisboa, Lisboa, Portugal

###License
This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License. To view a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/.
