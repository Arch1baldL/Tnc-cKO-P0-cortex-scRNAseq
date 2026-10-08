#!/bin/bash
#SBATCH --job-name=KO-2_Raw
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cellranger count \
  --id=KO-2 \
  --fastqs="$HOME/KO-2_raw" \
  --sample=KO-2 \
  --transcriptome="$HOME/refdata-gex-GRCm39-2024-A" \
  --create-bam=true \
  --force-cells=12000