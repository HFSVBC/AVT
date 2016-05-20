# AVT - Allele Validation Tool (AVT) for whole genome Multilocus Sequence Typing (wgMLST) schemas
AVT (Allele Validation Tool), is a tool capable of validate novel alleles found in de novo assemblers by allele calling algorithms.

The tool needs four / five inputs. A txt file with all the inf and references (example: example/AlleleCallresults.txt). A second txt file with the nodes and it's refrences (example: example/ContigLocation.txt). A third file with the sequence (example: example/10_S10_L001.fasta). And the read(s) (example: example/10_S10_L001_R1_001/ and example/10_S10_L001_R2_001/). The output will be a txt file with the alleles that have not passed the validation.

AVT's website: http://hfsvbc.github.io/AVT/

**Example**
```
python AVT.py AlleleCallresults.txt ContigLocation.txt 10_S10_L001.fasta 10_S10_L001_R1_001.fastq -R2 10_S10_L001_R2_001.fastq -d True
```

**Example of Output File**
```
#CHROM	POS	ID	REF	ALT	ERROR_PERCENTAGE(ALT_is_Wrong)
gi_407930685_ref_NC_018706.1_:c230209-229088.fasta	1430	.	G	T	3.15573117368e-21
gi_407930685_ref_NC_018706.1_:c230209-229088.fasta	1208	.	G	A	7.92683831666e-16
gi_384141246_ref_NC_017171.1_:2598773-2599849.fasta	58	.	C	A	3.97283016734e-16
gi_407930685_ref_NC_018706.1_:c230209-229088.fasta	1211	.	T	A	1.25631920858e-17
gi_384141246_ref_NC_017171.1_:2598773-2599849.fasta	355	.	G	A	3.15573117368e-21
gi_384141246_ref_NC_017171.1_:2598773-2599849.fasta	49	.	T	A	1.25631920858e-15
gi_565636615_ref_NC_023028.1_:3720022-3722859.fasta	2120	.	T	C	3.15573117368e-21
gi_384129960_ref_NC_017162.1_:c35008-33635.fasta	697	.	A	C	3.15573117368e-21
gi_407930685_ref_NC_018706.1_:c230209-229088.fasta	1409	.	C	T	3.15573117368e-21
gi_384141246_ref_NC_017171.1_:95696-96739.fasta	1297	.	G	A	3.15573117368e-21
```

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
This work is licensed under the GNU GPL V3.0 license, for more information check license.txt.
