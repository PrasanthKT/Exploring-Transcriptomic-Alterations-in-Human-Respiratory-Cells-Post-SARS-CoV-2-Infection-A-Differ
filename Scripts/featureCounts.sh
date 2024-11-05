#!/bin/bash

#SBATCH --job-name=featureCounts
#SBATCH --output=featureCounts_output_%j.txt
#SBATCH --error=featureCounts_error_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=06:00:00
#SBATCH --mem=100GB
#SBATCH --partition=general
#SBATCH -A r00750

# Load the module  miniconda
module load miniconda

# Activate the conda environement that contains the tool subread for the feature counts
conda activate snakemake_env

# Path to the annotation file
ANNOTATION="/N/u/pthuthi/Quartz/Assignment_2/Homo_sapiens.GRCh38.109.gtf"

# Path to the base directory that contains the SRR directories
BASE_DIR="/N/u/pthuthi/Quartz/Assignment_2"

# Find all sorted BAM files within the `mapping_results` directories
BAM_FILES=$(find "$BASE_DIR" -type f -path "*/mapping_results/*" -name "*_aligned_sorted.bam")

# Check if BAM files were found
if [ -z "$BAM_FILES" ]; then
  echo "No BAM files found. Exiting..."
  exit 1
fi

# Run featureCounts on all found BAM files and save it to .csv format
featureCounts -a "$ANNOTATION" -o "$BASE_DIR/feature_counts.csv" $BAM_FILES

echo "featureCounts completed. Results saved to feature_counts.csv"
