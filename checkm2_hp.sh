#!/bin/bash
# =============================================================================
# CheckM2 Quality Assessment of H. pylori Genomes
#
# This script runs CheckM2 on all downloaded H. pylori genome assemblies
# to assess completeness and contamination before downstream analysis.
#
# We process genomes in batches of 200 because loading the full DIAMOND
# database into memory requires ~16GB RAM. Running smaller batches keeps
# memory usage within a standard 8GB workstation.
#
# Input:  hp_fasta_files_all/   — folder of .fna genome assemblies
# Output: checkm2_results/      — per-batch quality reports
#         quality_report_all.tsv — combined report for all genomes
#
# Tools required:
#   CheckM2 v1.0.1 (conda: checkm2)
#   Python 3 with pandas
#
# Author: Hasnain Jadoon
# Institution: Dalian University of Technology / BGI Genomics, Shenzhen
# Study: Multi-species Comparative Genomics of Gastric Cancer-Associated Bacteria
# =============================================================================

set -euo pipefail

GENOME_DIR="hp_fasta_files_all"
RESULTS_DIR="checkm2_results"
COMBINED_TSV="quality_report_all.tsv"
BATCH_SIZE=200

# Activate environment
conda activate checkm2

# Check input directory exists
if [ ! -d "$GENOME_DIR" ]; then
    echo "Error: Genome directory '$GENOME_DIR' not found."
    echo "Run download_hp_genomes.sh first to download genome assemblies."
    exit 1
fi

total=$(ls "$GENOME_DIR"/*.fna 2>/dev/null | wc -l)
echo "Found $total genome assemblies in $GENOME_DIR"
echo "Processing in batches of $BATCH_SIZE to manage memory usage"
echo ""

# Split genomes into batch folders
echo "Setting up batch directories..."
ls "$GENOME_DIR"/*.fna | split -l "$BATCH_SIZE" - genome_batch_list_

batch_num=0
for batch_list in genome_batch_list_*; do
    batch_num=$((batch_num + 1))
    batch_dir="checkm2_batch_${batch_num}"
    result_dir="${RESULTS_DIR}/batch_${batch_num}"

    # Skip if this batch already has results
    if [ -f "${result_dir}/quality_report.tsv" ]; then
        echo "Batch $batch_num already complete, skipping"
        rm -f "$batch_list"
        continue
    fi

    # Create symlinks to genome files for this batch
    mkdir -p "$batch_dir"
    while read fna_path; do
        ln -sf "$(realpath "$fna_path")" "$batch_dir/"
    done < "$batch_list"

    genome_count=$(ls "$batch_dir"/*.fna | wc -l)
    echo "Running CheckM2 on batch $batch_num ($genome_count genomes)..."

    checkm2 predict \
        --input "$batch_dir/" \
        --output-directory "$result_dir" \
        --extension fna \
        --threads 4 \
        --lowmem \
        --force

    echo "  Batch $batch_num complete"

    # Clean up batch directory and list
    rm -rf "$batch_dir"
    rm -f "$batch_list"
done

echo ""
echo "All batches complete. Combining results..."

# Combine all batch reports into one file
python3 << 'PYEOF'
import pandas as pd
import glob

reports = sorted(glob.glob("checkm2_results/batch_*/quality_report.tsv"))
print(f"Found {len(reports)} batch reports")

combined = pd.concat([pd.read_csv(f, sep='\t') for f in reports], ignore_index=True)
combined.to_csv("quality_report_all.tsv", sep='\t', index=False)

total = len(combined)
hq = combined[(combined['Completeness'] >= 90) & (combined['Contamination'] <= 5)]
low = combined[combined['Completeness'] < 90]
cont = combined[combined['Contamination'] > 5]

print(f"\n=== Quality Assessment Summary ===")
print(f"Total genomes assessed                     : {total}")
print(f"High quality (>=90% complete, <=5% contam) : {len(hq)} ({len(hq)/total*100:.1f}%)")
print(f"Failed completeness (<90%)                 : {len(low)}")
print(f"Failed contamination (>5%)                 : {len(cont)}")
print(f"\nResults saved to: quality_report_all.tsv")
PYEOF

echo ""
echo "Done. Quality report saved to: $COMBINED_TSV"
echo "Use this report to filter genomes before pangenome analysis."
