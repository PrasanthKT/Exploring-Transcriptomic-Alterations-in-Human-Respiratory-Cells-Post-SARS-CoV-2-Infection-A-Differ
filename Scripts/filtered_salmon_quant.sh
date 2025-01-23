#!/bin/bash
#SBATCH -J salmon_quant_filtered
#SBATCH -p general
#SBATCH -o salmon_quant_filtered_%j.out
#SBATCH -e salmon_quant_filtered_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=8:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

module load miniconda

conda activate rna-seq

BASE_DIR="/N/scratch/pthuthi/RNA-Seq_check"
INDEX_DIR="/N/scratch/pthuthi/RNA_Seq/indexing/salmon_index"
OUTPUT_DIR="/N/scratch/pthuthi/RNA_Seq/quantification/salmon_filtered"

mkdir -p "${OUTPUT_DIR}"

RNA_SEQ_RANGE=$(seq -f "SRR12375%03g" 92 106)

for SRR in $RNA_SEQ_RANGE; do
    FILTER_DIR="${BASE_DIR}/${SRR}/prinseq_filtered"
    
    if [ -d "$FILTER_DIR" ]; then
        FILTERED_FASTQ="${FILTER_DIR}/${SRR}_filtered.fastq"
        if [ -f "$FILTERED_FASTQ" ]; then
            SAMPLE_OUTPUT="${OUTPUT_DIR}/${SRR}"
            mkdir -p "$SAMPLE_OUTPUT"
            echo "Running Salmon quantification for filtered file: $FILTERED_FASTQ..."
            salmon quant -i "$INDEX_DIR" \
                         -l A \
                         -r "$FILTERED_FASTQ" \
                         -p 8 \
                         -o "$SAMPLE_OUTPUT"
            if [ $? -eq 0 ]; then
                echo "Salmon quantification completed for $SRR."
            else
                echo "Salmon quantification failed for $SRR."
            fi
        else
            echo "Filtered FASTQ file not found for $SRR: $FILTERED_FASTQ"
        fi
    else
        echo "Filtered directory not found for $SRR: $FILTER_DIR"
    fi
done

echo "Salmon quantification for filtered files completed."
