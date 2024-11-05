#!/bin/bash

# Adapter sequence for trimming # According to the fastqc report it says it has Universal illumina adapter content so using that sequence to remove the adapters in the reads.
ADAPTER="AGATCGGAAGAG"

# Base directory where all SRR directories are located
BASE_DIR="/N/u/pthuthi/Quartz/Assignment_2"

# Loop through each SRR directory and run Cutadapt
for DIR in "$BASE_DIR"/SRR*; do
  if [ -d "$DIR" ]; then
    # Identify the .fastq file in the directory
    for FASTQ_FILE in "$DIR"/*.fastq; do
      if [ -f "$FASTQ_FILE" ]; then
        # Define output file path
        OUTPUT_FILE="${FASTQ_FILE%.fastq}_trimmed.fastq"
        
        echo "Running Cutadapt on $FASTQ_FILE..."
        
        # Run Cutadapt for adapter trimming
        cutadapt -a "$ADAPTER" -o "$OUTPUT_FILE" "$FASTQ_FILE"
        
        echo "Trimming completed for $FASTQ_FILE, output saved to $OUTPUT_FILE"
      else
        echo "No .fastq files found in $DIR"
      fi
    done
  fi
done

echo "Adapter trimming completed for all files."
