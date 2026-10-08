#!/bin/bash
#SBATCH --job-name=Ctrl-2_Raw
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cellranger count \
  --id=Ctrl-2 \
  --fastqs="$HOME/Ctrl-2_raw" \
  --sample=Ctrl-2 \
  --transcriptome="$HOME/refdata-gex-GRCm39-2024-A" \
  --create-bam=true \
  --force-cells=10000