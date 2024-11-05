#!/bin/bash

# Base directory that contains all the fastq files directories
BASE_DIR="/N/u/pthuthi/Quartz/Assignment_2"

# Loop over each directory in the base directory
for DIR in "$BASE_DIR"/SRR*; do
  if [ -d "$DIR" ]; then
    # Create a 'fastqc' directory in the current directory to store FastQC results
    FASTQC_DIR="$DIR/fastqc"
    mkdir -p "$FASTQC_DIR"

    # Loop through each .fastq file in the directory
    for FASTQ_FILE in "$DIR"/*.fastq; do
      if [ -f "$FASTQ_FILE" ]; then
        echo "Running FastQC on $FASTQ_FILE..."
        
        # Run FastQC and save the results in the fastqc directory
        fastqc -o "$FASTQC_DIR" "$FASTQ_FILE"
        
        echo "Finished FastQC on $FASTQ_FILE"
      else
        echo "No .fastq files found in $DIR"
      fi
    done
  fi
done

echo "FastQC analysis completed for all .fastq files."

