#!/bin/bash
#SBATCH -J hisat2_mapping
#SBATCH -p general
#SBATCH -o hisat2_mapping_%j.txt
#SBATCH -e hisat2_mapping_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --time=48:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

module load miniconda

conda activate wgcna

BASE_DIR="/N/scratch/pthuthi/wgcna"
REFERENCE_DIR="/N/scratch/pthuthi/wgcna/reference"
GENOME_FA="$REFERENCE_DIR/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
FILTERED_SUFFIX="prinseq_filtered"

echo "Building HISAT2 index..."
hisat2-build -p 4 "$GENOME_FA" "$REFERENCE_DIR/GRCh38_index"
if [ $? -ne 0 ]; then
    echo "Error: HISAT2 index building failed."
    exit 1
fi
echo "HISAT2 index built successfully."

for SRR_ID in $(seq -f "SRR12375%03g" 92 106); do
    SRR_DIR="$BASE_DIR/$SRR_ID"
    if [ -d "$SRR_DIR" ]; then
        FASTQ_DIR="$SRR_DIR/$FILTERED_SUFFIX"
        for FASTQ_FILE in "$FASTQ_DIR"/*.fastq; do
            if [ -f "$FASTQ_FILE" ]; then
                SAMPLE_NAME=$(basename "$FASTQ_FILE" .fastq)
                OUTPUT_SAM="$SRR_DIR/${SAMPLE_NAME}_aligned.sam"
                OUTPUT_BAM="$SRR_DIR/${SAMPLE_NAME}_aligned_sorted.bam"

                echo "Mapping $FASTQ_FILE to reference genome..."
                hisat2 -p 4 -x "$REFERENCE_DIR/GRCh38_index" -U "$FASTQ_FILE" -S "$OUTPUT_SAM"

                if [ $? -eq 0 ]; then
                    echo "Converting SAM to sorted BAM for $SAMPLE_NAME..."
                    samtools view -bS "$OUTPUT_SAM" | samtools sort -o "$OUTPUT_BAM"
                    samtools index "$OUTPUT_BAM"

                    # Cleanup intermediate SAM file to save space
                    rm "$OUTPUT_SAM"
                    echo "Mapping and sorting completed for $SAMPLE_NAME."
                else
                    echo "Error: Mapping failed for $SAMPLE_NAME. Check the log for details."
                fi
            else
                echo "No FASTQ files found in $FASTQ_DIR."
            fi
        done
    else
        echo "Directory $SRR_DIR not found for $SRR_ID."
    fi
done

echo "HISAT2 mapping completed for all RNA-Seq files in the specified range."
