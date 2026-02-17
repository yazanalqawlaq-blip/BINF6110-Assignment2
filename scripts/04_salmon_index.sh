#!/usr/bin/env bash
#SBATCH --job-name=salmon_index
#SBATCH --account=def-cottenie
#SBATCH --time=00:30:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm/04_salmon_index_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL
set -euo pipefail
cd /scratch/yazanalq/BINF6110-Assignment2

# Build Salmon index with genome as decoy (prevents genomic reads from
# mismapping to the transcriptome — see Srivastava et al. 2020)

module load salmon/1.10.2

REF_DIR=reference
INDEX_DIR=salmon_index

# Prepare the decoy: extract chromosome names from the genome
grep "^>" ${REF_DIR}/GCF_000146045.2_R64_genomic.fna | cut -d " " -f 1 | sed 's/>//' > ${REF_DIR}/decoys.txt

# Concatenate transcriptome + genome (transcriptome must come first)
cat ${REF_DIR}/GCF_000146045.2_R64_rna.fna ${REF_DIR}/GCF_000146045.2_R64_genomic.fna > ${REF_DIR}/gentrome.fna

# Build the index
salmon index \
    -t ${REF_DIR}/gentrome.fna \
    -d ${REF_DIR}/decoys.txt \
    -i ${INDEX_DIR} \
    -p 8

echo "Salmon index built — $(date)"
ls -lh ${INDEX_DIR}/
