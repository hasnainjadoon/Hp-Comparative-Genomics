#!/bin/bash
# ═══════════════════════════════════════════════════════════════
# Download All H. pylori Genomes from NCBI via CLI
# Total unique genomes: 9,385 (deduplicated GCA+GCF)
# ═══════════════════════════════════════════════════════════════

# ── STEP 0: Install NCBI datasets CLI ─────────────────────────
conda install -c conda-forge ncbi-datasets-cli -y

# Or download latest binary directly (if conda network is blocked)
mkdir -p ~/bin
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/datasets -O ~/bin/datasets
wget https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/LATEST/linux-amd64/dataformat -O ~/bin/dataformat
chmod +x ~/bin/datasets ~/bin/dataformat
echo 'export PATH="$HOME/bin:$PATH"' >> ~/.bashrc
source ~/.bashrc

# ── STEP 1: Download Metadata for ALL H. pylori Genomes ───────
# Taxon ID 210 = Helicobacter pylori
# Includes all assembly levels (complete, chromosome, scaffold, contig)
datasets summary genome taxon 210 \
    --assembly-level complete,chromosome,scaffold,contig \
    --as-json-lines > hp_all_genomes_metadata.jsonl

# ── STEP 2: Convert JSON to Clean TSV Table ───────────────────
dataformat tsv genome \
    --inputfile hp_all_genomes_metadata.jsonl \
    --fields accession,organism-name,organism-tax-id,source_database,assminfo-level,assminfo-refseq-category,assminfo-release-date,assminfo-biosample-accession,assmstats-total-sequence-len,assmstats-number-of-contigs,assmstats-contig-n50,assmstats-contig-l50,assmstats-gc-percent,checkm-completeness,checkm-contamination,assminfo-biosample-strain,assminfo-biosample-host,assminfo-biosample-host-disease,assminfo-biosample-geo-loc-name,assminfo-biosample-isolation-source \
    > hp_all_genomes_metadata.tsv

echo "Total genome records: $(wc -l < hp_all_genomes_metadata.tsv) (including header)"

# ── STEP 3: Deduplicate GCA/GCF Paired Records ────────────────
# Many genomes appear twice (once as GCA GenBank + once as GCF RefSeq)
# This keeps one record per BioSample, preferring RefSeq (GCF) when available
python3 << 'EOF'
import pandas as pd

df = pd.read_csv('hp_all_genomes_metadata.tsv', sep='\t')
print(f"Before deduplication: {len(df)} records")
print(f"  GCA (GenBank): {df['Assembly Accession'].str.startswith('GCA_').sum()}")
print(f"  GCF (RefSeq) : {df['Assembly Accession'].str.startswith('GCF_').sum()}")

# Sort so GCF comes first, then deduplicate by BioSample
df['is_refseq'] = df['Assembly Accession'].str.startswith('GCF_')
df_sorted = df.sort_values('is_refseq', ascending=False)
df_dedup = df_sorted.drop_duplicates(subset='Assembly BioSample Accession', keep='first')
df_dedup = df_dedup.drop(columns=['is_refseq'])

print(f"After deduplication : {len(df_dedup)} unique genomes")
df_dedup.to_csv('hp_all_genomes_deduplicated.tsv', sep='\t', index=False)
print("Saved: hp_all_genomes_deduplicated.tsv")
EOF

# ── STEP 4: Check Gastric Cancer Strains ──────────────────────
echo ""
echo "=== Host Disease Breakdown ==="
cut -f18 hp_all_genomes_deduplicated.tsv | sort | uniq -c | sort -rn | head -30

echo ""
echo "=== Gastric Cancer Strains ==="
cut -f18 hp_all_genomes_deduplicated.tsv | \
    grep -i "cancer\|carcinoma\|adenocarcinoma\|malignan\|neoplas\|tumor\|tumour" | \
    sort | uniq -c | sort -rn

# Extract GC strain list
awk -F'\t' '
NR==1 {print; next}
tolower($18) ~ /cancer|carcinoma|adenocarcinoma|malignan|neoplas|tumor|tumour/ && $18 !~ /non cancer/ {print}
' hp_all_genomes_deduplicated.tsv > hp_gc_strains_final.tsv

echo ""
echo "GC strains saved: $(( $(wc -l < hp_gc_strains_final.tsv) - 1 )) genomes → hp_gc_strains_final.tsv"

# ── STEP 5: Extract Accession List for Download ───────────────
cut -f1 hp_all_genomes_deduplicated.tsv | tail -n +2 > hp_all_9385_accessions.txt
echo "Accession list saved: $(wc -l < hp_all_9385_accessions.txt) accessions"

# ── STEP 6: Split Into Batches of 500 ─────────────────────────
split -l 500 hp_all_9385_accessions.txt hp_batch_
echo "Batches created: $(ls hp_batch_?? | wc -l)"

# ── STEP 7: Download Using Dehydrate/Rehydrate Method ─────────
# NOTE: Standard streaming downloads fail in China due to network issues
# Dehydrate downloads a tiny manifest first, then rehydrates files individually
# This is far more resilient to connection drops

mkdir -p hp_fasta_files_all

for batch in hp_batch_??; do
    dehydrated_zip="${batch}_dehydrated.zip"
    dehydrated_dir="${batch}_dehydrated"

    # Skip if already fully rehydrated
    if [ -d "$dehydrated_dir" ] && \
       [ "$(find "$dehydrated_dir" -name '*.fna' | wc -l)" -gt 0 ]; then
        echo "✔ $batch already done ($(find "$dehydrated_dir" -name '*.fna' | wc -l) genomes)"
        continue
    fi

    echo "► Downloading manifest: $batch..."
    datasets download genome accession --inputfile "$batch" \
        --include genome \
        --dehydrated \
        --filename "$dehydrated_zip"

    if [ -f "$dehydrated_zip" ]; then
        mkdir -p "$dehydrated_dir"
        unzip -q "$dehydrated_zip" -d "$dehydrated_dir"

        echo "  Rehydrating $batch (downloading genome files)..."
        cd "$dehydrated_dir"
        datasets rehydrate --directory .
        cd ..

        count=$(find "$dehydrated_dir" -name "*.fna" | wc -l)
        echo "  ✅ $batch: $count genomes downloaded"
    else
        echo "  ❌ $batch: manifest download failed"
        echo "$batch" >> failed_batches.txt
    fi
done

# ── STEP 8: Collect All FASTA Files Into One Folder ───────────
echo ""
echo "Collecting all .fna files..."
for dir in hp_batch_*_dehydrated; do
    find "$dir" -name "*.fna" -exec cp {} hp_fasta_files_all/ \;
done

total=$(ls hp_fasta_files_all/*.fna | wc -l)
echo "✅ Total genomes collected: $total"
echo "Location: hp_fasta_files_all/"

# ── STEP 9: Optional Cleanup ──────────────────────────────────
# Only run after confirming all 9385 .fna files are in hp_fasta_files_all/
# rm -rf hp_batch_*_dehydrated hp_batch_*_dehydrated.zip

# ── FINAL SUMMARY ─────────────────────────────────────────────
echo ""
echo "════════════════════════════════════════"
echo "SUMMARY"
echo "════════════════════════════════════════"
echo "Total H. pylori genomes (raw)    : $(( $(wc -l < hp_all_genomes_metadata.tsv) - 1 ))"
echo "After deduplication (unique)     : $(( $(wc -l < hp_all_genomes_deduplicated.tsv) - 1 ))"
echo "Gastric cancer-source strains    : $(( $(wc -l < hp_gc_strains_final.tsv) - 1 ))"
echo "FASTA files downloaded           : $(ls hp_fasta_files_all/*.fna 2>/dev/null | wc -l)"
echo "Output folder                    : hp_fasta_files_all/"
