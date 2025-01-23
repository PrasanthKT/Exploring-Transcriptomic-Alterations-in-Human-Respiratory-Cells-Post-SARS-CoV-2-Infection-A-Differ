#!/bin/bash
#SBATCH -J salmon_quant
#SBATCH -p general
#SBATCH -o salmon_quant_%j.out
#SBATCH -e salmon_quant_%j.err
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
INDEX_DIR="${BASE_DIR}/indexing/salmon_index"
OUTPUT_DIR="${BASE_DIR}/quantification/salmon"

mkdir -p "${OUTPUT_DIR}"

for SAMPLE_DIR in ${BASE_DIR}/SRR*; do
    SAMPLE_ID=$(basename "${SAMPLE_DIR}")
    FASTQ_FILE="${SAMPLE_DIR}/${SAMPLE_ID}.fastq"
    
    if [[ -f "$FASTQ_FILE" ]]; then
        echo "Running Salmon quantification for ${SAMPLE_ID}..."
        salmon quant -i "${INDEX_DIR}" \
                     -l A \
                     -r "${FASTQ_FILE}" \
                     -p 8 \
                     -o "${OUTPUT_DIR}/${SAMPLE_ID}"
        echo "Salmon quantification completed for ${SAMPLE_ID}."
    else
        echo "FASTQ file missing for ${SAMPLE_ID}, skipping..."
    fi
done

echo "Salmon quantification completed for all samples."
