#!/bin/bash
#SBATCH --job-name=Ctrl-3_Raw
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cellranger count \
  --id=Ctrl-3 \
  --fastqs="$HOME/Ctrl-3_raw" \
  --sample=Ctrl-3 \
  --transcriptome="$HOME/refdata-gex-GRCm39-2024-A" \
  --create-bam=true \
  --force-cells=12000