#!/bin/bash

# Base directory that contains all the Raw SRA files downloaded using SRA_Tool kit
BASE_DIR="/N/u/pthuthi/Quartz/Assignment_2"

# Loop over each directory in the base directory
for DIR in "$BASE_DIR"/SRR*; do
  if [ -d "$DIR" ]; then
    # Loop through each .sra file in the directory
    for SRA_FILE in "$DIR"/*.sra; do
      if [ -f "$SRA_FILE" ]; then
        echo "Processing $SRA_FILE..."
        
        # Run fastq-dump for single-end files
        fastq-dump --outdir "$DIR" --skip-technical --readids --dumpbase --clip "$SRA_FILE"
        
        echo "Finished processing $SRA_FILE"
      else
        echo "No .sra files found in $DIR"
      fi
    done
  fi
done

echo "All .sra files have been processed."

