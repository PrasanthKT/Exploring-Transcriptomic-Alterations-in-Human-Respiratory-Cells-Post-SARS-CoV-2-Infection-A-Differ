#!/bin/bash
#SBATCH -J Mapping_the_reads_        
#SBATCH -p general                       
#SBATCH -o mapping_output_%j.txt     
#SBATCH -e mapping_error%j.err      
#SBATCH --mail-type=ALL                  
#SBATCH --mail-user=pthuthi@iu.edu       
#SBATCH --nodes=1                        
#SBATCH --ntasks-per-node=1              
#SBATCH --time=12:00:00                 
#SBATCH --mem=200GB                       
#SBATCH -A r00750

# Load the module load miniconda
module load miniconda

#Activate the conda environment that cntains the hisat2 aligner
conda activate rna_env

#Note: Index the reference genome before.

# Directory that contains Reference genome to map the reads
REFERENCE="/N/u/pthuthi/Quartz/Assignment_2/Inputs/Homo_sapiens_index"

# Base directory containing input files
BASE_DIR="/N/u/pthuthi/Quartz/Assignment_2"

# Loop through each directory with adapater_trimmed FASTQ files
for DIR in "$BASE_DIR"/SRR*; do
  if [ -d "$DIR" ]; then
    # Identify the trimmed FASTQ file
    for FASTQ_FILE in "$DIR"/*_trimmed.fastq; do
      if [ -f "$FASTQ_FILE" ]; then
        # Create a directory for the mapping results
        RESULT_DIR="$DIR/mapping_results"
        mkdir -p "$RESULT_DIR"
        
        # Define the output SAM file path
        SAM_FILE="$RESULT_DIR/$(basename "$FASTQ_FILE" _trimmed.fastq)_aligned.sam"
        
        echo "Mapping $FASTQ_FILE to $REFERENCE..."
        
        # Run HISAT2 for mapping
        hisat2 -p 4 -x "$REFERENCE" -U "$FASTQ_FILE" -S "$SAM_FILE"
        
        echo "Mapping completed for $FASTQ_FILE. Results saved in $SAM_FILE"
      else
        echo "No trimmed FASTQ file found in $DIR"
      fi
    done
  fi
done

echo "All mappings completed." 
