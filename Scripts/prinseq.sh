#!/bin/bash
#SBATCH -J preprocess_rnaseq
#SBATCH -p general
#SBATCH -o preprocess_rnaseq_%j.txt
#SBATCH -e preprocess_rnaseq_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=20:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

module load miniconda

conda activate wgcna

BASE_DIR="/N/scratch/pthuthi/RNA-Seq_check"

RNA_SEQ_RANGE=$(seq -f "SRR12375%03g" 92 106)

for SRR in $RNA_SEQ_RANGE; do
    DIR="$BASE_DIR/$SRR"
    if [ -d "$DIR" ]; then
        
        FILTER_DIR="$DIR/prinseq_filtered"
        mkdir -p "$FILTER_DIR"

        for FASTQ_FILE in "$DIR"/*.fastq; do
            if [ -f "$FASTQ_FILE" ]; then
                echo "Running PRINSEQ on RNA-Seq file: $FASTQ_FILE..."
                prinseq-lite.pl -fastq "$FASTQ_FILE" \
                                -out_good "$FILTER_DIR/${SRR}_filtered" \
                                -out_bad "$FILTER_DIR/${SRR}_discarded" \
                                -derep 1 \
                                -lc_method dust \
                                -lc_threshold 7
                if [ $? -eq 0 ]; then
                    echo "PRINSEQ completed for $FASTQ_FILE."
                else
                    echo "PRINSEQ failed for $FASTQ_FILE."
                fi
            fi
        done
    else
        echo "Directory $DIR not found."
    fi
done

echo "RNA-Seq preprocessing completed."

