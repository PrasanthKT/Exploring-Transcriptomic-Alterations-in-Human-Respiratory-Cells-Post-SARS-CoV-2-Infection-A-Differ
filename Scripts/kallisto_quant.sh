#!/bin/bash
#SBATCH -J kallisto_quant
#SBATCH -p general
#SBATCH -o kallisto_quant_%j.out
#SBATCH -e kallisto_quant_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

module load miniconda

conda activate rna-seq

BASE_DIR="/N/scratch/pthuthi/RNA_Seq"
INDEX_DIR="${BASE_DIR}/indexing/kallisto_index/kallisto.idx"
OUTPUT_DIR="${BASE_DIR}/quantification/kallisto"


mkdir -p "${OUTPUT_DIR}"

for SAMPLE_DIR in ${BASE_DIR}/SRR*; do
    SAMPLE_ID=$(basename "${SAMPLE_DIR}")
    FASTQ_FILE="${SAMPLE_DIR}/${SAMPLE_ID}.fastq"
    
    if [[ -f "$FASTQ_FILE" ]]; then
        echo "Running Kallisto quantification for ${SAMPLE_ID}..."
        kallisto quant -i "${INDEX_DIR}" \
                       -o "${OUTPUT_DIR}/${SAMPLE_ID}" \
                       --single -l 200 -s 20 \
                       -t 8 \
                       "${FASTQ_FILE}"
        echo "Kallisto quantification completed for ${SAMPLE_ID}."
    else
        echo "FASTQ file missing for ${SAMPLE_ID}, skipping..."
    fi
done

echo "Kallisto quantification completed for all samples."
