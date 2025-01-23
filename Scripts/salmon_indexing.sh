#!/bin/bash
#SBATCH -J salmon_indexing
#SBATCH -p general
#SBATCH -o salmon_indexing_%j.out
#SBATCH -e salmon_indexing_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

module load miniconda

conda activate rna-seq

BASE_DIR="/N/scratch/pthuthi/RNA_Seq/reference"
TRANSCRIPTOME="${BASE_DIR}/Homo_sapiens.GRCh38.cdna.all.fa"
OUTPUT_DIR="/N/scratch/pthuthi/RNA_Seq/indexing/salmon_index"

mkdir -p "${OUTPUT_DIR}"

echo "Starting Salmon Indexing..."
salmon index -t "${TRANSCRIPTOME}" \
             -i "${OUTPUT_DIR}" \
             -k 31 \
             --keepDuplicates

if [ $? -eq 0 ]; then
    echo "Salmon Indexing Completed Successfully."
else
    echo "Salmon Indexing Failed."
    exit 1
fi
