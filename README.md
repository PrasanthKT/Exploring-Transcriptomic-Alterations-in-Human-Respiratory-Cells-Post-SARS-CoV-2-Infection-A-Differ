# Differential Gene Expression analysis in PM2.5-Exposed Bronchial Epithelial Cells

### Objective
To uncover the molecular mechanisms underlying PM2.5 exposure in bronchial epithelial cells by analyzing transcriptomic data. This study employs differential gene expression analysis to identify key genes, regulatory pathways, and biological processes disrupted by environmental pollutants.

### Tools and Dependencies
1. SRA-Toolkit
2. Fastqc
3. Prinseq
4. Trim Galore
5. HISAT2
6. Sam-tools
7. Subread
8. DESeq2
9. Cluster Profiler
10. R and Bio conductor Packages
11. GO

### Installation

Create a conda environment.
```
conda create -n rna-seq
conda activate rna-seq
```
Install the dependencies using bioconda
```
conda install -c bioconda sra-tools fastqc trim-galore prinseq hisat2 samtools subread
```
For R and Bioconductor packages, install them directly in R:
```
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "ComplexHeatmap"))
```

### Input files 
1. SRA IDs: RNA-Seq Data: SRP275647 Bisulfite - Seq Data: SRP275645
2. Reference Genome: Homo_sapiens.GRCh38.dna.primary_assembly.fa
3. Annotation File: Homo_sapiens.GRCh38.109.gtf

1. Download Raw data using SRA Toolkit
```
prefetch SRP275647
prefetch SRP275645
```
2. Human Reference Genome 
```
wget ftp://ftp.ensembl.org/pub/release-109/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```
3. Humman Annotation File
```
wget ftp://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/Homo_sapiens.GRCh38.109.gtf.gz
gunzip Homo_sapiens.GRCh38.109.gtf.gz
```
### Getting Started
1.	Set up the Conda environment and install all dependencies.
2.	Download the raw SRA files, reference genome, and annotation file using the commands above.
3.	Navigate to the Scripts folder and execute the preprocessing scripts for RNA-Seq and Bisulfite-Seq data.
4.	Use the output files from the analysis for further processing and analysis in R.

### Workflow
The detailed workflow for the analysis is documented in the report, which is available in the repository as ```Report.pdf```. Shell scripts for each step of the workflow are organized under the ```Scripts``` folder, while R scripts and their corresponding input files are located in the ```DESeq2```folder within the repository.

### Results
A comprehensive interpretation of the results is provided in the ```Report.pdf``` report. Key findings are summarized, with visualizations and pathway insights derived from the integrative analysis of transcriptomic and epigenomic data.
