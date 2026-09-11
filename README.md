# Hp-Comparative-Genomics

Bioinformatics pipeline for downloading and analyzing *Helicobacter pylori* genomes from NCBI for gastric cancer pathogenesis research.

## Overview
- Downloads all H. pylori genomes (taxon ID: 210) from NCBI using datasets CLI
- Deduplicates GCA/GCF paired records → 9,385 unique genomes
- Extracts metadata including host disease annotation
- Identifies 218 confirmed gastric cancer-source strains
- Downloads FASTA files using dehydrate/rehydrate method (resilient to network drops)

## Key Results
| Metric | Count |
|--------|-------|
| Total NCBI records | 13,066 |
| Unique genomes (deduplicated) | 9,385 |
| Gastric cancer-source strains | 182 |

## Requirements
- NCBI datasets CLI (v18.32+)
- Python 3 with pandas

## Usage
```bash
chmod +x download_hp_genomes.sh
./download_hp_genomes.sh 2>&1 | tee download_hp.log
```

## Author
Hasnain Jadoon
Master's researcher, Dalian University of Technology / BGI Genomics, Shenzhen
Published: Gut 2026;75(Suppl 2):A262-A264
DOI: https://doi.org/10.1136/gutjnl-2026-IDDF.166
