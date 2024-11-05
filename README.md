## miRNA_Analysis

### Objective 
The primary goal of this project is to analyze transcriptomic changes in human respiratory cells following SARS-CoV-2 infection. This involves identifying differentially expressed genes (DEGs) and performing gene ontology (GO) and pathway enrichment analyses to understand the biological processes and pathways affected by the virus at different time points.

### Data Source
Data for this project is available through the National Center for Biotechnology Information (NCBI) under the accession BioProjectID: PRJNA901149.

### Tools/Dependences
1. Conda
2. Bioconda
3. SRA-Toolkit
4. Fastqc
5. Cutadapt
6. Hisat2
7. Sam-tools
8. Subread
9. DESQ2

### Installation
1. Create a conda environment for the project
```bash
conda create -n name_of_the_env
```
2. Install the tools using the following commands in the created environments.
```bash
conda install -c bioconda sra-tools fastqc cutadapt hisat2 samtools subread
```
3. For DESQ2 in R
```R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
```
### Download the data and preprocess
```bash
prefetch SRP407642
```
Note: These are single ended reads
```bash
fastq-dump --split-files
```
### Reference Genome and Annotation file
1. Download the reference genome file 
```bash
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
2. Download the gene annotation file 
```bash
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz 
``` 
### Analysis
The next scripts of the script can be acessed under ```Scripts``` folder present in the repository. Make sure you change your paths of the directories and the environment names according to your system paths and modify the slurm details to submit the job according to your slurm account details.
```fastq-dump.sh``` To convert the .sra files to fastq files
```fastqc.sh```  Quality Control for the fastq files to know the quality of the fastq files.
```cutadapt.sh``` After analysing the fastq report it shows that it contains the adapter content, it has illuminated universal adapter content ```AGATCGGAAGAG``` so used the cutadapt was used to remove the adapters.
```mapping.sh``` Mapping the reads to the refernce genome
```sam_to_bam.sh``` Convert the .same file format to .bam format for the further analysis.
```featureCounts.sh``` Using subread to do the feature counts using the gtf file.
### Differential gene expression analysis.
Detailed R script for the differential gene expression analysis is presnt in the repository

### Output files
The output files for each step are used for the analysis of the next steps. Results from all the steps are summarized and the detailed report ```Report.pdf``` in the repository.

