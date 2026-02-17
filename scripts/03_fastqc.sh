#!/usr/bin/env bash
#SBATCH --job-name=fastqc
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --output=slurm/03_fastqc_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL
set -euo pipefail
cd /scratch/yazanalq/BINF6110-Assignment2

# Run FastQC on all raw reads
module load fastqc/0.12.1

fastqc raw_data/*.fastq.gz \
    --outdir fastqc_results \
    --threads 4

echo "FastQC complete — $(date)"
ls -lh fastqc_results/
