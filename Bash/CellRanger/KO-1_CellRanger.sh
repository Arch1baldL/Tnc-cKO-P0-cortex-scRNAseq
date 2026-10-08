#!/bin/bash
#SBATCH --job-name=KO-1_Raw
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cellranger count \
  --id=KO-1 \
  --fastqs="$HOME/KO-1_raw" \
  --sample=KO-1 \
  --transcriptome="$HOME/refdata-gex-GRCm39-2024-A" \
  --create-bam=true \
  --force-cells=10000