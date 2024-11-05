#!/bin/bash

# Base directory that contains all the SRR directories
BASE_DIR="/N/u/pthuthi/Quartz/Assignment_2"

# Loop through each SRR directory and process SAM files in `mapping_results`
for DIR in "$BASE_DIR"/SRR*/mapping_results; do
  if [ -d "$DIR" ]; then
    # Process each SAM file in the mapping_results directory
    for SAM_FILE in "$DIR"/*.sam; do
      if [ -f "$SAM_FILE" ]; then
        # Define output BAM file path
        BAM_FILE="${SAM_FILE%.sam}.bam"

        echo "Converting $SAM_FILE to BAM format..."
        
        # Convert SAM to BAM
        samtools view -bS "$SAM_FILE" > "$BAM_FILE"
        
        echo "Conversion completed for $SAM_FILE, output saved as $BAM_FILE"
      else
        echo "No SAM files found in $DIR"
      fi
    done
  fi
done

echo "All SAM to BAM conversions completed."

