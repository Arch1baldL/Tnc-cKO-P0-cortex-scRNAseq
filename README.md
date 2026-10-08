# Tnc-cKO-P0-cortex-scRNAseq
Code and analysis workflows for single-cell RNA-seq of the P0 mouse dorsal cortex following *Tnc* conditional knockout.

<br>

## Overview
This repository contains the analysis code accompanying the study:

**A Single-Cell Transcriptomic Dataset of the Neonatal Mouse Dorsal Cortex Following Tenascin-C Deletion**

The dataset comprises single-cell RNA-sequencing data from the dorsal cortex of postnatal day 0 (P0) mice, including 3 control (Ctrl) and 3 *Tnc* conditional knockout (KO) samples.

<br>

## Repository Structure
```text
Tnc-cKO-P0-cortex-scRNAseq/
├── Bash/
│   ├── CellRanger/
│   │   ├── Ctrl-1_CellRanger.sh
│   │   ├── Ctrl-2_CellRanger.sh
│   │   ├── Ctrl-3_CellRanger.sh
│   │   ├── KO-1_CellRanger.sh
│   │   ├── KO-2_CellRanger.sh
│   │   └── KO-3_CellRanger.sh
│   └── velocyto/
│       ├── Ctrl-1_velocyto.sh
│       ├── Ctrl-2_velocyto.sh
│       ├── Ctrl-3_velocyto.sh
│       ├── KO-1_velocyto.sh
│       ├── KO-2_velocyto.sh
│       └── KO-3_velocyto.sh
├── R/
│   ├── 01_Integration.R
│   └── 02_Annotation.R
├── Python/
│   ├── 01_RNA-velocity.py
│   └── 02_Velocity-QC.py
├── LICENSE
└── README.md
```

<br>

Scripts within each directory are numbered according to the recommended analysis workflow.

<br>

## Data availability
Raw sequencing data and the processed Seurat object are available from Science Data Bank:

**[A Single-Cell Transcriptomic Dataset of the Neonatal Mouse Dorsal Cortex Following Tenascin-C Deletion](https://doi.org/10.57760/sciencedb.013su)**

<br>

## Reproducibility
The scripts in this repository reflect the analysis workflow used for the associated study.
Paths to large input datasets and computing environments may need to be adjusted before execution on another system. Parameters affecting the reported analyses are retained in the corresponding scripts.

<br>

## License
The code in this repository is released under the MIT License.
The sequencing and processed datasets are distributed separately through Science Data Bank and are not covered by the software license of this repository.
