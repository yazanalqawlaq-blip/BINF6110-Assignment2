#!/usr/bin/env bash
#SBATCH --job-name=salmon_quant
#SBATCH --account=def-cottenie
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --output=slurm/05_salmon_quant_%j.out
#SBATCH --mail-user=yazanahmadalqawlaq@gmail.com
#SBATCH --mail-type=END,FAIL
set -euo pipefail
cd /scratch/yazanalq/BINF6110-Assignment2

# Quantify all 9 samples with Salmon
# Single-end reads require --fldMean and --fldSD (estimated fragment length)
# Using typical Illumina RNA-seq values: mean=180, sd=20

module load salmon/1.10.2

INDEX_DIR=salmon_index
RAW_DIR=raw_data
OUT_DIR=salmon_output

SAMPLES=(
    SRR10551665
    SRR10551664
    SRR10551663
    SRR10551662
    SRR10551661
    SRR10551660
    SRR10551659
    SRR10551658
    SRR10551657
)

for SRR in "${SAMPLES[@]}"; do
    echo "Quantifying: ${SRR} — $(date)"

    salmon quant \
        -i ${INDEX_DIR} \
        -l A \
        -r ${RAW_DIR}/${SRR}.fastq.gz \
        -p 8 \
        --fldMean 180 \
        --fldSD 20 \
        --validateMappings \
        --gcBias \
        -o ${OUT_DIR}/${SRR}

    echo "${SRR} done."
done

echo ""
echo "Quantification summary:"
for SRR in "${SAMPLES[@]}"; do
    MAPPED=$(grep "num_mapped" ${OUT_DIR}/${SRR}/aux_info/meta_info.json | tr -dc '0-9.')
    echo "${SRR}: ${MAPPED} mapped reads"
done
