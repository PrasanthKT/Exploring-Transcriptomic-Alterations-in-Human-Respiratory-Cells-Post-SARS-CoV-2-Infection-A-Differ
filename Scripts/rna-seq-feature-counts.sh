#!/bin/bash
#SBATCH -J featureCounts
#SBATCH -p general
#SBATCH -o featureCounts_%j.txt
#SBATCH -e featureCounts_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=pthuthi@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --time=4:00:00
#SBATCH --mem=20GB
#SBATCH -A r00750

module load miniconda

conda activate wgcna

BASE_DIR="/N/scratch/pthuthi/wgcna"
FEATURECOUNTS_RESULTS="$BASE_DIR/Feature_Counts_All"
REFERENCE_GTF="$BASE_DIR/reference/Homo_sapiens.GRCh38.109.gtf"

mkdir -p "$FEATURECOUNTS_RESULTS"

BAM_FILES=""
for SRR_ID in SRR12375092 SRR12375093 SRR12375094 SRR12375095 SRR12375096 SRR12375097 SRR12375098 SRR12375099 SRR12375100 SRR12375101 SRR12375102 SRR12375103 SRR12375104 SRR12375105 SRR12375106; do
    BAM_FILE="$BASE_DIR/$SRR_ID/${SRR_ID}_filtered_aligned_sorted.bam"
    if [ -f "$BAM_FILE" ]; then
        BAM_FILES="$BAM_FILES $BAM_FILE"
    else
        echo "Warning: BAM file $BAM_FILE not found. Skipping..."
    fi
done

if [ -n "$BAM_FILES" ]; then
    echo "Running featureCounts..."
    featureCounts -T 4 -a "$REFERENCE_GTF" -o "$FEATURECOUNTS_RESULTS/feature_counts.txt" $BAM_FILES


    if [ -f "$FEATURECOUNTS_RESULTS/feature_counts.txt" ]; then
        echo "Converting feature_counts.txt to feature_counts.csv..."
        awk -v OFS="," 'BEGIN {FS="\t"} NR==2 {print $0} NR>2 {print $0}' "$FEATURECOUNTS_RESULTS/feature_counts.txt" > "$FEATURECOUNTS_RESULTS/feature_counts.csv"
        echo "Feature counts saved to $FEATURECOUNTS_RESULTS/feature_counts.csv"
    else
        echo "Error: feature_counts.txt not found. Conversion to CSV failed."
    fi

    echo "Feature counting completed successfully."
else
    echo "No BAM files found. Exiting..."
    exit 1
fi
