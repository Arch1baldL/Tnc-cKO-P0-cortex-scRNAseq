#!/bin/bash
#SBATCH --job-name=Ctrl-1_Raw
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32

cellranger count \
  --id=Ctrl-1 \
  --fastqs="$HOME/Ctrl-1_raw" \
  --sample=Ctrl-1 \
  --transcriptome="$HOME/refdata-gex-GRCm39-2024-A" \
  --create-bam=true