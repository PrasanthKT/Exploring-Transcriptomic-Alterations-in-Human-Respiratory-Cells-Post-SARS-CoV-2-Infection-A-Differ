#!/bin/bash
#SBATCH -J kallisto_indexing
#SBATCH -p general
#SBATCH -o kallisto_indexing_%j.out
#SBATCH -e kallisto_indexing_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=10:00:00
#SBATCH --mem=20GB
#SBATCH -A r00750

module load miniconda

conda activate rna-seq

BASE_DIR="/N/scratch/pthuthi/RNA_Seq/reference"
TRANSCRIPTOME="${BASE_DIR}/Homo_sapiens.GRCh38.cdna.all.fa"
OUTPUT_DIR="/N/scratch/pthuthi/RNA_Seq/indexing/kallisto_index"

mkdir -p "${OUTPUT_DIR}"

echo "Starting Kallisto Indexing..."
kallisto index -i "${OUTPUT_DIR}/kallisto.idx" "${TRANSCRIPTOME}"

if [ $? -eq 0 ]; then
    echo "Kallisto Indexing Completed Successfully."
else
    echo "Kallisto Indexing Failed."
    exit 1
fi
