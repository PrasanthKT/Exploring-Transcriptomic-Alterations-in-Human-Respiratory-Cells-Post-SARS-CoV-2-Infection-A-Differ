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

