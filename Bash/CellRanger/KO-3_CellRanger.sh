#!/bin/bash
#SBATCH --job-name=KO-3_Raw
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cellranger count \
  --id=KO-3 \
  --fastqs="$HOME/KO-3_raw" \
  --sample=KO-3 \
  --transcriptome="$HOME/refdata-gex-GRCm39-2024-A" \
  --create-bam=true