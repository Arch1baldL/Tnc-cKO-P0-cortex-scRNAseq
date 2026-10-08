#!/bin/bash
#SBATCH --job-name=velocyto_KO-3
#SBATCH --partition=premium
#SBATCH --nodes=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G

set -euo pipefail

SAMPLE_ID="KO-3"
SAMPLE_DIR="${HOME}/${SAMPLE_ID}"
GTF_FILE="${HOME}/refdata-gex-GRCm39-2024-A/genes/genes.gtf"

export MAMBA_ROOT_PREFIX="${HOME}/micromamba"
source "${MAMBA_ROOT_PREFIX}/etc/profile.d/micromamba.sh"
micromamba activate velocyto_env

THREADS="${SLURM_CPUS_PER_TASK:-10}"

echo "Job Start Time: $(date)"
echo "Running on node: $(hostname)"
echo "Using velocyto:"
which velocyto
velocyto --version || true

echo "Starting velocyto run10x..."

velocyto run10x \
    -@ "${THREADS}" \
    --samtools-memory 4000 \
    "${SAMPLE_DIR}" \
    "${GTF_FILE}"

echo "Job End Time: $(date)"