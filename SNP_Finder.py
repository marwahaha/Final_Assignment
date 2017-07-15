import subprocess
import sys
import re
import glob


def file_exist(file):
    ''' Search the specified directory for a file'''
    files = glob.glob(file + '*')
    return len(files) > 0


def download_ref_genome(url):
    ''' Downloads a file from the provided url using the wget function'''
    # launch the subprocess wget in the comand line to retreive a file from the provided url
    p = subprocess.Popen(['wget',  url], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)


def unzip_ref_genome(path):
    ''' Unzip the file provided to path '''
    # launch the subprocess gunzip
    p = subprocess.Popen(['gunzip', path], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)


def is_fasta(filename, maxline=100):
    '''Determine if the file is a valid FASTA file with nucleotide sequences '''
    pattern = re.compile('^[ACGTWRKYSMBDHVN?-]+$')  # regular expression for nucleotide sequence
    nrecords = 0
    with open(filename, 'rU') as handle:
        for ln, line in enumerate(handle):
            if line.startswith('>'):
                nrecords += 1
                continue
            if not pattern.findall(line.upper()) and len(line.strip()) > 0:
                # non-header lines should contain valid nucleotide sequence
                return False
            if ln > maxline:
                break
    if nrecords == 0:
        return False
    return True


def make_reference(path):
    ''' Build a indexed reference file using Bowtie2 build function '''
    # pass this comand to the console bowtie2-build GCA_001624185.1_129S1_SvImJ_v1_genomic.fna --threads 8 musRef
    p = subprocess.Popen(['bowtie2-build', path, '--threads 8', 'musRef'], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)


def align_sam(reference, fastq1, fastq2):
    ''' Align the reference genome created by bowtie2 with paired fastq files. This function is meant for paired end reads'''
    #pass this comand to the console bowtie2 -x musRef -1 SRR3168572_1.fastq  -2 SRR3168572_2.fastq  -S mus_aligned --thread 8
    mus_aligned = 'mus_aligned'
    p = subprocess.Popen(['bowtie2', '-x', reference, '-1', fastq1, '-2', fastq2, '-S', mus_aligned, '--thread 8'], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)
    return mus_aligned


def sam2bam(path):
    ''' Convert a SAM file format into a BAM file format'''
    # use the comand samtools view -ub -h mus_aligned -o m_alig.bam -@8
    mus_alig = 'mus_alig.bam'
    subprocess.check_call(['samtools', 'view', '-ub', '-h', path, '-o', mus_alig, '-@8'], stdout=subprocess.PIPE)
    return mus_alig


def sort_bam(bam_file_path):
    ''' Sorts the BAM file based on the leftmost variable, usually the chromosome name'''
    # uses the command samtools sort -T aln.sorted -o mus_alig.sorted.bam -@8 mus_alig.bam
    sorted_bam = bam_file_path.replace(".bam", ".sorted.bam")
    p = subprocess.Popen(['samtools', 'sort', '-T',  'alng.sorted.tmp', '-o', sorted_bam, '-@8', bam_file_path], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)

    return sorted_bam


def index_sorted_bam(sorted_bam_path):
    ''' Creates an index for a sorted BAM file '''
    # samtools index mus_alig.sorted.bam
    p = subprocess.Popen(['samtools', 'index', sorted_bam_path], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)


def region_filter(sorted_indexed_bam):
    ''' Isolates the chromosomes and regions of interest to avoid working with the whole genome'''
    # uses the command samtools view -ub -h -@8 mus_alig.sorted.bam -o chrome_filtered.sorted.bam CM003937.1 CM003939.1 CM003941.1
    chrome_filtered = 'chrom_filter.sorted.bam'
    p = subprocess.Popen(['samtools', 'view', '-ub', '-h', '-@8', sorted_indexed_bam, '-o', chrome_filtered, 'CM003937.1', 'CM003939.1', 'CM003941.1'], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)
    return chrome_filtered


def make_VCF(refrence,chrom_filtered):
    ''' Create a vcf file using the samtools mpileup command and the bcftools call command '''
    # uses the command samtools mpileup -uf GCA_001624185.1_129S1_SvImJ_v1_genomic.fna chrom_filtered.sorted.bam | bcftools call -mv > var.raw.vcf
    p = subprocess.Popen(['samtools', 'mpileup', '-uf', refrence, chrom_filtered, '|', 'bcftools', 'call', '-mv', '>', 'chromosomes.vcf'], stdout=subprocess.PIPE)
    for line in p.stdout:
        ln = line.decode('ascii')
        print(ln)
    return vcf_file


def parse_VCF(path):
    ''' Extract the Chromosome,Position, Reference Allele, Alternate Allele, and Genotype from a VCF file.
     To change the chromosome and position needed the chromosomes dictionary needs to be modified with the required information '''
    # Create a dictionary with the chromosome name as depicted in the vcf file and the chromosome position of interest
    chromosomes = {"CM003937.1": ["56534778", "55427885"], "CM003939.1": ["143810388"], "CM003941.1": ["28970903"]}

    # open the input and output files
    handle = open(path, "rU")
    outfile = open("SNPFinder_results.csv", "w")

    # Create the header for the output
    header = ",".join(["Chromosome", "Position", "Reference Allele", "Alternate Allele", "Genotype", "\n"])
    outfile.write(header)

    for line in handle:
        # Skip the comment lines in the vcf file
        if line.startswith("##"):
            continue
        else:
            # turn the line into a tuple
            tuples = line.split("\t")
            chromo = tuples[0]
            pos = tuples[1]
            # Isolate the chromosome name and the position if they match the dictionary values
            if chromo in chromosomes.keys() and pos in chromosomes[chromo]:
                # Isolate the other variables of interest from the rest of the line
                RefAllele = tuples[3]
                AltAllele = tuples[4]
                genotypes = tuples[9]
                Genotype = genotypes[0:3]

                # Bring all the variables of interest together in a csv file
                FinalLine = ",".join([chromo, pos, RefAllele, AltAllele, Genotype, '\n'])
                outfile.write(FinalLine)
    handle.close()
    outfile.close()


def error_exit(msg):
    ''' Prints error message specified and exits the function '''
    print(msg)
    sys.exit(1)


BASE_DIR = "."
MUS_REFERENCE = "musRef*"
GENOME_FASTA_GZ_URL = "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/624/185/GCA_001624185.1_129S1_SvImJ_v1/GCA_001624185.1_129S1_SvImJ_v1_genomic.fna.gz"
GENOME_FASTA_GZ = 'GCA_001624185.1_129S1_SvImJ_v1_genomic.fna.gz'
GENOME_FASTA = 'GCA_001624185.1_129S1_SvImJ_v1_genomic.fna'

GENOME_FILE_ERROR = "Genome file is not FASTA"

# Step 1 : Get the reference file
print("Looking for reference files.....")
if not file_exist("/".join([BASE_DIR, MUS_REFERENCE])):
    # Check if fasta files exists
    if file_exist("/".join([BASE_DIR, GENOME_FASTA])):
        if is_fasta(GENOME_FASTA):
                make_reference(GENOME_FASTA)
        else:
            error_exit(GENOME_FILE_ERROR)
    else:
        # check if archive exits
        if file_exist("/".join([BASE_DIR, GENOME_FASTA_GZ])):
            print("Unzipping reference genome!")
            unzip_ref_genome(GENOME_FASTA_GZ)
            if is_fasta(GENOME_FASTA):
                make_reference(GENOME_FASTA)
            else:
                error_exit(GENOME_FILE_ERROR)
        else:
            # download the archive and unzip
            download_ref_genome(GENOME_FASTA_GZ_URL)
            print("Unzipping reference genome!")
            unzip_ref_genome(GENOME_FASTA_GZ)
            if is_fasta(GENOME_FASTA):
                make_reference(GENOME_FASTA)
            else:
                error_exit(GENOME_FILE_ERROR)


# Step 2 :
print("Aligning fastq file to reference genome ...... This will also take a whhhhhhhhhile. Maybe go for a break while it does its thing")
sam_file = align_sam("musRef", "SRR3168572_1.fastq", "SRR3168572_2.fastq")

# Step 3 :
print("Converting the SAM file to a BAM file.")
bam_file = sam2bam(sam_file)


# Step 4 :
print("Sorting the Bam file.")
sorted_bam_file = sort_bam(bam_file)

# Step 5 :
print("Indexing the sorted file.")
index_sorted_bam(sorted_bam_file)

# Step 6 :
print("Getting rid of the chromosomes you dont need.")
Chromosomes = region_filter(sorted_bam_file)

# Step 7 :
print("Almost there! Comparing the aligned and indexed chromosomes of interest to the reference genome to find varients.")
vcf_file = make_VCF(GENOME_FASTA, Chromosomes)

# Step 8 :
print("Last step! Generating the parsed VCF file")
parse_VCF(vcf_file)

# Step 9 :
print("OMG WE'RE DONE! Check the directory for the SNPFinder_results.csv file")
