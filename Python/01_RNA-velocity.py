# Setup and Configuration
from pathlib import Path

import anndata
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv

ROOT = Path(__file__).resolve().parent
DATASET = ROOT / "data"
INPUT = ROOT / "input"
OUT = ROOT / "output"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
SAMPLES = ["Ctrl1", "Ctrl2", "Ctrl3", "KO1", "KO2", "KO3"]
N_TOP_GENES = 3000
N_NEIGHBORS = 15
CELLTYPE_ORDER = [
    "RG", "NIPC", "GIPC", "OPC", "ASC", "EN", "IN", "CR",
    "MG", "Endo", "Fibro", "Epend", "RBC", "lowQC",
]
CELLTYPE_COLORS = {
    "RG": "#6A51A3", "GIPC": "#D7301F", "OPC": "#F5B70A",
    "ASC": "#FCAE6B", "NIPC": "#315BA8", "EN": "#7FB3D5",
    "IN": "#00B8A9", "CR": "#C44E52", "MG": "#009E73",
    "Endo": "#7F9A3A", "Epend": "#B5B35C", "Fibro": "#8E6C8A",
    "RBC": "#8C3B3B", "lowQC": "#9E9E9E",
}

mpl.rcParams.update({
    "font.family": "Liberation Sans",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
})

OUT.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

# Load input data
meta = pd.read_csv(INPUT / "TNC_umap_metadata.csv", index_col="barcode")
harmony = pd.read_csv(INPUT / "TNC_harmony_coords.csv", index_col=0)
pca = pd.read_csv(INPUT / "TNC_pca_coords.csv", index_col=0)
# Match cells and combine samples
sample_objects = []
match_rows = []
for sample in SAMPLES:
    loom_path = DATASET / sample / f"{sample}.loom"
    loom = sc.read_loom(loom_path, sparse=True, cleanup=False)
    loom.var_names_make_unique()
    names = loom.obs_names.str.replace(":", "_", regex=False).str.replace("x$", "-1", regex=True)
    prefix = f"{sample}_"
    loom.obs_names = pd.Index([name if name.startswith(prefix) else prefix + name for name in names])
    expected = pd.Index(meta.index[meta["orig.ident"].astype(str).eq(sample)])
    matched = expected.intersection(loom.obs_names)
    match_rows.append({
        "sample": sample,
        "tnc_cells": len(expected),
        "loom_cells": loom.n_obs,
        "matched_cells": len(matched),
        "tnc_match_rate": len(matched) / len(expected),
    })
    if len(matched) != len(expected):
        raise ValueError(f"{sample}: matched {len(matched)}/{len(expected)} TNC cells.")
    sample_objects.append(loom[matched].copy())

pd.DataFrame(match_rows).to_csv(TABLES / "barcode_match.csv", index=False)
adata = anndata.concat(sample_objects, axis=0, join="inner", index_unique=None)
# Prepare AnnData object
adata = adata[meta.index].copy()
adata.obs = meta.loc[adata.obs_names].copy()
adata.obs["cell_type"] = pd.Categorical(
    adata.obs["cell_type"].astype(str), categories=CELLTYPE_ORDER, ordered=True
)
adata.uns["cell_type_colors"] = [CELLTYPE_COLORS[label] for label in CELLTYPE_ORDER]
adata.obs["condition"] = np.where(adata.obs["orig.ident"].astype(str).str.startswith("Ctrl"), "Ctrl", "KO")
adata.obsm["X_harmony"] = harmony.loc[adata.obs_names].to_numpy(dtype=np.float32)
adata.obsm["X_pca"] = pca.loc[adata.obs_names].to_numpy(dtype=np.float32)
adata.obsm["X_umap"] = meta.loc[adata.obs_names, ["UMAP_1", "UMAP_2"]].to_numpy(dtype=np.float32)
adata.X = (adata.layers["spliced"] + adata.layers["unspliced"]).tocsr()

# Filter and select genes
scv.pp.filter_and_normalize(adata, min_shared_counts=20)
score = np.asarray(adata.layers["spliced"].mean(axis=0)).ravel()
score += np.asarray(adata.layers["unspliced"].mean(axis=0)).ravel()
keep = np.argsort(score)[-N_TOP_GENES:]
adata = adata[:, np.sort(keep)].copy()

# Calculate RNA velocity
for key in ["neighbors", "pca"]:
    adata.uns.pop(key, None)
for key in ["distances", "connectivities"]:
    adata.obsp.pop(key, None)
sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X_harmony", random_state=123)
scv.pp.moments(adata, n_pcs=None, n_neighbors=N_NEIGHBORS)
scv.tl.velocity(adata, mode="stochastic")
scv.tl.velocity_graph(adata, n_jobs=4, backend="threading", show_progress_bar=False)
scv.tl.velocity_confidence(adata)
scv.tl.velocity_pseudotime(adata)
scv.tl.velocity_embedding(adata, basis="umap")

# Save velocity object
adata.write_h5ad(OUT / "velocity.h5ad", compression="gzip")


# Save gene-level QC
gene_qc = adata.var[["velocity_genes", "velocity_gamma", "velocity_r2"]].copy()
gene_qc.index.name = "gene"
gene_qc = gene_qc.reset_index()
gene_qc.to_csv(TABLES / "gene_qc.csv", index=False)

# Plot velocity grid
scv.set_figure_params("scvelo", dpi=120, dpi_save=300)
grid_fig, grid_ax = plt.subplots(figsize=(8, 8))
scv.pl.velocity_embedding_grid(
    adata, basis="umap", color="cell_type",
    palette=[CELLTYPE_COLORS[label] for label in CELLTYPE_ORDER], density=1.0, arrow_size=1.2,
    title="TNC RNA velocity", legend_loc="right margin",
    show=False, ax=grid_ax,
)
grid_ax.set_aspect("equal", adjustable="box")
grid_fig.savefig(FIGURES / "umap_grid.png", dpi=400)
grid_fig.savefig(FIGURES / "umap_grid.pdf")
plt.close("all")

# Plot velocity stream
scv.pl.velocity_embedding_stream(
    adata, basis="umap", color="cell_type",
    palette=[CELLTYPE_COLORS[label] for label in CELLTYPE_ORDER], density=1.0,
    title="TNC RNA velocity", legend_loc="right margin", show=False,
)
plt.savefig(FIGURES / "umap_stream.png", dpi=400, bbox_inches="tight")
plt.savefig(FIGURES / "umap_stream.pdf", bbox_inches="tight")
plt.close("all")
