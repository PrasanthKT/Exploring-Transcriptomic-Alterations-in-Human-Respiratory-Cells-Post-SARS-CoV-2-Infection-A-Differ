# Assessing the Impact of Different Tools on Transcriptomic Analysis Methods

### Objective
This project compares transcriptomic workflows with varierty of tools at each step under various computational constraints, using PM_2.5 exposed bronchial epithelial cell data. The goal is to see whether workflow selection and tool slection significantly influences the final results or not, compared to the original study.

### Tools and Dependencies
1. SRA-Toolkit
2. Fastqc
3. Prinseq
4. Salmon
5. Kallisto
6. DESeq2
7. EdgeR
8. R and Bio conductor Packages

### Installation

Create a conda environment.
```
conda create -n rna-seq
conda activate rna-seq
```
Install the dependencies using bioconda
```
conda install -c bioconda sra-tools fastqc prinseq salmon kallisto
```
For R and Bioconductor packages, install them directly in R:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "edgeR"))
```

### Input files 
1. SRA IDs: RNA-Seq Data: SRP275647
2. Reference Transciptome: Homo_sapiens.GRCh38.cdna.all.fa
3. Annotation File: Homo_sapiens.GRCh38.109.gtf

1. Download Raw data using SRA Toolkit
```
prefetch SRP275647
```
2. Human Reference Transcriptome 
```
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
3. Humman Annotation File
```
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz
```
### Getting Started
1.	Set up the Conda environment and install all dependencies.
2.	Download the raw SRA files, reference transcriptome, and annotation file using the commands above.
3.	Navigate to the Scripts folder and execute scripts for RNA-Seq within the conda environment.
4.	Use the output files from the analysis for further processing and analysis in R.

### Workflow
The detailed workflow for the analysis is documented in the report, which is available in the repository as ```Report.pdf```. Shell scripts for each step of the workflow are organized under the ```Scripts``` folder, while R scripts and their corresponding input files are located in the ```R_Analysis```folder within the repository.

### Results
A comprehensive interpretation of the results is provided in the ```Report.pdf``` report. Key findings are summarized in the rport.
