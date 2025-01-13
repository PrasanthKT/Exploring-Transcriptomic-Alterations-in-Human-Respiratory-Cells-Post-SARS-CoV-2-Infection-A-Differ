#!/bin/bash
#SBATCH -J fastqc_post_prinseq
#SBATCH -p general
#SBATCH -o fastqc_post_prinseq_%j.txt
#SBATCH -e fastqc_post_prinseq_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

module load miniconda

conda activate wgcna

BASE_DIR="/N/scratch/pthuthi/wgcna"

RNA_SEQ_RANGE=$(seq -f "SRR12375%03g" 92 106)

for SRR in $RNA_SEQ_RANGE; do
    DIR="$BASE_DIR/$SRR"
    if [ -d "$DIR" ]; then
    
        FILTERED_DIR="$DIR/prinseq_filtered"
        FASTQC_DIR="$DIR/processed_fastqc"
 
        mkdir -p "$FASTQC_DIR"

        for FASTQ_FILE in "$FILTERED_DIR"/*.fastq; do
            if [ -f "$FASTQ_FILE" ]; then
                echo "Running FastQC on file: $FASTQ_FILE..."
                fastqc -o "$FASTQC_DIR" "$FASTQ_FILE"
            else
                echo "No FASTQ files found in $FILTERED_DIR for $SRR."
            fi
        done
    else
        echo "Directory $DIR not found for $SRR."
    fi
done

echo "FastQC analysis completed for all preprocessed RNA-Seq files."
