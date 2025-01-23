#!/bin/bash
#SBATCH -J fastqc
#SBATCH -p general
#SBATCH -o fastqc_%j.txt
#SBATCH -e fastqc_%j.err
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

        REPORT_DIR="$DIR/fastqc_reports"
        mkdir -p "$REPORT_DIR"
        
        for FASTQ_FILE in "$DIR"/*.fastq; do
            if [ -f "$FASTQ_FILE" ]; then
                echo "Running FastQC on $FASTQ_FILE..."
                fastqc "$FASTQ_FILE" --outdir="$REPORT_DIR"
                
                if [ $? -eq 0 ]; then
                    echo "FastQC completed for $FASTQ_FILE. Results saved in $REPORT_DIR."
                else
                    echo "Error running FastQC on $FASTQ_FILE."
                fi
            else
                echo "No FASTQ files found in $DIR."
            fi
        done
    fi
done

echo "FastQC analysis completed for all FASTQ files."

