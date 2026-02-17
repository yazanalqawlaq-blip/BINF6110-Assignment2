#!/usr/bin/env bash
#SBATCH --job-name=download
#SBATCH --account=def-cottenie
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=slurm/01_download_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL
set -euo pipefail
cd /scratch/yazanalq/BINF6110-Assignment2/raw_data

# Download raw RNA-seq reads from NCBI SRA
# Mardanov et al. (2020) - PRJNA592304
# S. cerevisiae flor yeast, single-end Illumina
# 9 samples: 3 biofilm stages x 3 replicates

module load sra-toolkit/3.0.9

SAMPLES=(
    SRR10551665  # Early biofilm (IL20)
    SRR10551664  # Early biofilm (IL21)
    SRR10551663  # Early biofilm (IL22)
    SRR10551662  # Thin biofilm (IL23)
    SRR10551661  # Thin biofilm (IL24)
    SRR10551660  # Thin biofilm (IL25)
    SRR10551659  # Mature biofilm (IL29)
    SRR10551658  # Mature biofilm (IL30)
    SRR10551657  # Mature biofilm (IL31)
)

for SRR in "${SAMPLES[@]}"; do
    echo "Processing: ${SRR} — $(date)"

    if [ -f "${SRR}.fastq.gz" ]; then
        echo "${SRR} already exists, skipping."
        continue
    fi

    prefetch ${SRR} --max-size 5G
    fasterq-dump ${SRR}/${SRR}.sra --split-3 -e 1
    gzip ${SRR}.fastq
    rm -rf ${SRR}/
    echo "${SRR} done."
done

echo ""
echo "Download summary:"
ls -lh *.fastq.gz
echo "Expected: 9 | Found: $(ls *.fastq.gz | wc -l)"
