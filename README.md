# Somatic Variant Analysis in Mouse Glioma Tumor Samples Using Whole Genome Sequencing Data and Genome Analysis Tool Kit (GATK)

### Objective
This project involves performing a complete NGS pipeline for variant calling, starting from downloading the raw sequencing data to final variant annotation. Additionally, functional enrichment analysis is conducted to understand the biological significance of the variants identified in two types of mouse glioma tumor samples.

### Tools and Dependencies
1. SRA-Toolkit
2. Fastqc
3. Trim Galore
4. BWA
5. Sam-tools
6. GATK
7. Vcf-tools
8. Picard
9. ANNOVAR
10. ClueGO

### Installation

Create a conda environment.
```
conda create -n gatk
```
Install the dependencies using the following command.
```
conda install -c bioconda sra-tools
conda install -c bioconda fastqc
conda install -c bioconda trim-galore
conda install -c bioconda bwa
conda install -c bioconda samtools
conda install -c bioconda gatk
conda install -c bioconda vcftools
conda install -c bioconda picard
```
### Input files 
1. SRA IDs: The sequencing data for three samples were downloaded using the SRA Toolkit.
2. Reference Genome: Mus_musculus.GRCm39
3. Annotation File: Mus_musculus.GRCm39.109.gff3
