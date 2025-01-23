#!/bin/bash
#SBATCH -J fasterq-dump
#SBATCH -p general
#SBATCH -o fasterq-dump_%j.txt
#SBATCH -e fasterq-dump_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=2:00:00
#SBATCH --mem=10GB
#SBATCH -A r00750

module load miniconda

conda activate rna-seq

BASE_DIR="/N/scratch/pthuthi/RNA_Seq"

for DIR in $BASE_DIR/SRR*; do
    if [ -d "$DIR" ]; then
        
        cd "$DIR" || exit

        for SRA_FILE in *.sra; do
            if [ -f "$SRA_FILE" ]; then
                echo "Processing $SRA_FILE in $DIR..."

            
                fasterq-dump "$SRA_FILE" --threads 4 --outdir "$DIR"

                if [ $? -eq 0 ]; then
                    echo "$SRA_FILE successfully processed to FASTQ files."
                else
                    echo "Error processing $SRA_FILE. Check the log for details."
                fi
            else
                echo "No .sra files found in $DIR."
            fi
        done

        cd "$BASE_DIR" || exit
    fi
done

echo "fasterq-dump completed for all files."

