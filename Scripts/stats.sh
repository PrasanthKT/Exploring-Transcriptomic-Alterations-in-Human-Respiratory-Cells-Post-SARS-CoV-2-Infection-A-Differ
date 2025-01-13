#!/bin/bash

BASE_DIR="/N/scratch/pthuthi/wgcna"
RNA_SEQ_RANGE=$(seq -f "SRR12375%03g" 92 106)

for SRR_ID in $RNA_SEQ_RANGE; do
    BAM_FILE="$BASE_DIR/$SRR_ID/${SRR_ID}_filtered_aligned_sorted.bam"
    if [ -f "$BAM_FILE" ]; then
        echo "Mapping statistics for $SRR_ID:"
        samtools flagstat "$BAM_FILE"
        echo "----------------------------------"
    else
        echo "BAM file not found for $SRR_ID."
    fi
done
