# SNPFinder Pipeline #

**** Warning**: This script is designed to use 8 processors.

**** Warning**: This script starts from the very begining. And will take longer the first time it runs as it is creating the indexing files which will be used as it is run subsequently for other samples. There are checks in place to make sure it is not duplicating the most time consuming steps

## About the project ##
When viral infection strikes one of the organisms first lines of defense are the interferon induced restriction factors. These proteins interfere with several steps of the viral infection cycle including viral entry, replication of viral DNA, viral particle formation, viral budding and viral particle release. These genes have evolved to combat diverse viruses and inhibit their ability to replicate, propagate, and cause disease. The Barr lab studies disruption to the function of restriction factors that can arise due to a single base pair substitution in the DNA (*single nucleotide polymorphisms * **SNP**) resulting in non-synonoymous mutations or region deletions. These DNA mutations result in a loss of function for the restriction factors, which significantly influences disease progression. Our lab aims to develop a system whereby a profile of the mutations in the restriction factors of an individual can be correlated to disease progression. This is a pilot project to determine if we can identify the polymorphisms in a set of publicly available mouse genomes and scale it up to human datasets at a later time.

**Table 1**: Select SNP's in  129/SvlmJ mice complete with the SNP ID, chromosome position and Alleles

| Human Restriction Factor|Mouse Restriction Factor| Mouse SNP ID |Chromosome Location|Alleles|
| ------------------ |:-------------:| ----------:|---------------|-------------|
| HERC6              |HERC6| rs30287508  |Chr6:143810388| A/G|
| BST-2 / tetherin   |  BST-2   | rs50846085  |Chr8:71537427|C/T
|IFNβ | 	IFNβ1|rs28084066|Chr4:56534778|G/A|
|IFNα | IFNα 9|rs28083944|Chr4:55427885|G/T|

## Data required ##

1. Sample files to input into the script. This script takes paire end reads but can be modified to accept unpaired reads *see step 2*. As a proof of concept for the script I downloaded Illumina HiSeq 2000 paired end sequencing from the mouse tumor line ID8-G7 from the European Nucleotide Archive (http://www.ebi.ac.uk/ena/data/view/PRJNA310255)
	- SRR3168572_1.fastq
	- SRR3168572_2.fastq
	
2. Mouse Genome to create the refecrence file and align the sample sequences. If this is not provided the script will scrape the NCBI website and download the genome for Mus musculus	strain 129S1/SvImJ (https://www.ncbi.nlm.nih.gov/genome/genomes/52)

3. Bowtie 2 version 2.3.2 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)

4. samtools (Tools for alignments in the SAM format) Version: 1.5-2 (using htslib 1.5-5-gda5c0c7)

5. bcftools (Tools for data in the VCF/BCF formats) Version: 0.1.19


## Objectives

1. Create a reference genome using bowtie2
2. Map the sample files to the reference genome
3. Convert alligned sam file to the bianery bam format
4. Sort the bam file by read name
5. Index the sorted bam file
6. Isolate the chromosomes of interest
7. Create the VCF file using bcftools
8. Parse the vcf file to display the Gene, Chromosome Position,Reference Allele, Alternative Alleles, Genotyping (GT), Genotyping quality score (GQ)
9. Create the output .csv file

## Functions and their Descriptions

* `file_exist`: Search the specified directory for a file
	- input file name
	- input the directory (relative or absolute)
	- output returns true or false
* `download_ref_genome`: Downloads a file from the provided url using the wget function
	- input url
	- output the file specified by the url into the current directory
* `unzip_ref_genome`: Unzip the file provided to path
	- input file name if its in the current directory or the path(absolute or relative)
	- output unziped version of the input
* `is_fasta`: Determine if the file is a valid FASTA file with nucleotide sequences
	- input fastafile name if its in the current directory or file path (absolute or relative)
	- output return true or false
* `make_reference`: Build a indexed reference file using Bowtie2 build function 
	- input genome reference file in fasta format
	- output indexed reference genome in bowtie2 format (.bt2)
* `align_sam`: Align the reference genome created by bowtie2 with paired fastq files. This function is meant for paired end reads'
	- input bowtie2 formated reference files
	- input \*\_1.fastq sequence files
	- input \*\_2.fastq sequence files
	- output aligned SAM formated file
* `sam2bam`: Convert a SAM file format into a BAM file forma
	- input SAM formated file 
	- output bianery BAM formated file
* `sort_bam`: Sorts the BAM file based on the leftmost variable, usually the chromosome name
	- input aligned BAM file
	- output sorted BAM file (\*.sorted.bam)
* `index_sorted_bam`: Creates an index for a sorted BAM file
	- input sorted BAM file (\*.sorted.bam)
	- output index file in the current directory (.idx)
* `region_filter`: Isolates the chromosomes and regions of interest to avoid working with the whole genome
	- input aligned, sorted, indexed BAM file
	- output filtered BAM file containing only the chromosomes of interest
* `make_VCF`: Create a vcf file using the samtools mpileup command and the bcftools call command 
	- input filtered BAM file
	- output 
* `parse_VCF`: Extract the Chromosome,Position, Reference Allele, Alternate Allele, and Genotype from a VCF file.   To change the chromosome and position needed the chromosomes dictionary needs to be modified with the required information
	- input VCF formated file
	- output csv file
* `error_exit`: Prints error message specified and exits the function
	- input error message desired
    
    
## Steps and their input / output files ##
* Step 1: Search for the reference and the fasta files. If they are not there download and unzip.
	* Input files - bowtie2 indexed reference files OR fasta file of reference chromosome OR gziped fasta file of reference chromosome OR none
	* The script will run look in the current directory for the "musRef\*" bt2 index files. If the files are not there it will searched for the fasta file to create the reference from. If the fasta file is not there it will download and unzip the fasta file. If the fasta file or the gziped fasta file is in the current directory it will make the reference files.
	* reference file will be created using the bowtie2 command
```
bowtie2-build GCA_001624185.1_129S1_SvImJ_v1_genomic.fna --threads 8 musRef
```
	* GCA_001624185.1_129S1_SvImJ_v1_genomic.fna - is the fasta file
	* \-- threads 8 - uses 8 procesors to complete the task
	* musRef - is the output file names
	
* Step 2: Align the fastq sequence files with the indexed reference genome created in step 1.
```
bowtie2 -x musRef -1 SRR3168572_1.fastq  -2 SRR3168572_2.fastq  -S mus_aligned --thread 8
```
	* -x - The basename of the index for the reference genome.
	* musRef - the index genome in bt2 format
	* -1 - the first fastq file in a pair
	* -2 - th second fastq file in a pair
	* -S - output format is SAM
	* mus_aligned - is the output file name, this is a SAM file
	* \--threads 8 - is the number of procesors used

* Step 3: Convert the aligned file from step 2 into a binary BAM file for faster processing
```
samtools view -ub -h mus_aligned -o m_alig.bam -@8
```
	* -ub - the output file will be uncompressed bam 
	* -h - maintain header
	* -o - signifies output, output file name follows it
	* -@8 - using 8 procesors

* Step 4: Sort the BAM file from step 3 by the value in the lesfmost column 
```
samtools sort -T aln.sorted -o mus_alig.sorted.bam -@8 mus_alig.bam
```
	* -T - creates temporary file created so that everything is not loaded into memory
	* aln.sorted - temp file name
	* -o - the output file is BAM by default but need to specify the file name
	* mus\_alig.sorted.bam - output file name
	* -@8 - using 8 procesors
	* mus\_alig.bam - the input file name

* Step 5: Indexing the sorted file from step 4
```
samtools index mus\_alig.sorted.bam
```
	* input file is mus\_alig.sorted.bam this generates a new .idx file in the current directory with is used in the bacground but not directly at any step

* Step 6: Filter out the chromosomes you are not interested in and keep only the chromosomes that you are interested in
```
 samtools view -ub -h -@8 mus_alig.sorted.bam -o chrome_filtered.sorted.bam CM003937.1 CM003939.1 CM003941.1
```
	* -ub - the output file will be uncompressed bam 
	* -h - maintain header
	* -@8 - using 8 procesors
	* mus\_alig.sorted.bam - input file from step 5
	* -o - signifies output, output file name follows it
	* CM003937.1 CM003939.1 CM003941.1 - the names of the mouse chromosomes you need to keep. these chromosome names are based on the naming in the reference fasta file - can also be found in the website where the fasta file is obtaibed from
    
* Step 7: Align the procesed file with the original reference fasta file to determine the varience between the reference and the reads in the fastq files from step 2
```
samtools mpileup -uf GCA_001624185.1_129S1_SvImJ_v1_genomic.fna chrom_filtered.sorted.bam | bcftools call -mv > var.raw.vcf
```
	* -uf - the 'u' flag generates an uncompressed file (this was chosen for faster processing in the pipe to follow); the 'f' flag signals that the input file is in fasta format
	* GCA_001624185.1_129S1_SvImJ_v1_genomic.fna - the fasta file the reference genome was created from in step 1
	* chrom_filtered.sorted.bam - the input file from step 6 to align to the reference fasta file
	* | - pipes the output to bcftools
	* -mv - the 'm' flag is for multiallelic-caller this is an alternative model for multiallelic and rare-variant calling the 'v' flag outputs the varients only
	* \> - output is writen out to var.raw.vcf

* Step 8: Parses the VCF file from step 7 using the chromosome position
	* input file is the VCF file created in step 7 and the output file is a csv formated table. See Final Output file description below for details about the output information
	

## Final Output File description ##
* SNPFinder_results.csv
	- Gene - contains the abbreviated name of the gene of interest (ie. Hect and RLD containing protein 6 will be HERC6) 
    - Chromosome Position - position of the SNP on the chromosome
    - Reference Allele - the allele present in the reference genome
    - Alternative Alleles - the alternate allele present
    - Genotyping (GT) - homozygous for the reference allele (0/0), heterozygous (0/1), or homozygous for the aternative allele (1/1)
