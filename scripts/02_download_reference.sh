#!/usr/bin/env bash
# NOTE: Run on login node (not sbatch) — requires internet access
# Usage: bash scripts/02_download_reference.sh

set -euo pipefail
cd /scratch/yazanalq/BINF6110-Assignment2/reference

# Download S. cerevisiae R64 reference files from NCBI
BASE_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/146/045/GCF_000146045.2_R64"

wget ${BASE_URL}/GCF_000146045.2_R64_rna.fna.gz
wget ${BASE_URL}/GCF_000146045.2_R64_genomic.fna.gz
wget ${BASE_URL}/GCF_000146045.2_R64_genomic.gtf.gz

gunzip *.gz

echo "Reference files:"
ls -lh
