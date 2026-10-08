# Setup and Configuration
from pathlib import Path

import anndata
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scvelo as scv

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "output"
INPUT = OUT / "velocity.h5ad"
FIGURES = OUT / "figures"
TABLES = OUT / "tables"
FIGURES.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

mpl.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

# Plot functions
def save_figure(fig, stem):
    fig.savefig(FIGURES / f"{stem}.png", dpi=600)
    fig.savefig(FIGURES / f"{stem}.pdf")
    plt.close(fig)


def umap_panel(ax, xy, values, title, colorbar_label, cmap="viridis", vmin=None, vmax=None):
    finite = np.isfinite(xy).all(axis=1) & np.isfinite(values)
    points = ax.scatter(
        xy[finite, 0], xy[finite, 1], c=values[finite], s=0.55,
        cmap=cmap, vmin=vmin, vmax=vmax, linewidths=0, rasterized=True,
    )
    ax.set_title(title, pad=3)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_aspect("equal", adjustable="box")
    ax.tick_params(length=2, width=0.4, pad=1.5)
    for spine in ax.spines.values():
        spine.set_linewidth(0.4)
    colorbar = plt.colorbar(points, ax=ax, fraction=0.046, pad=0.025)
    colorbar.set_label(colorbar_label, fontsize=8)
    colorbar.ax.tick_params(labelsize=7, length=2, width=0.4)
    return points

# Load input data
adata = anndata.read_h5ad(INPUT)
# Fit dynamical model
scv.tl.recover_dynamics(
    adata,
    var_names="velocity_genes",
    n_jobs=4,
    backend="threading",
    show_progress_bar=False,
)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata, n_jobs=4, backend="threading", show_progress_bar=False)
scv.tl.velocity_confidence(adata)
scv.tl.latent_time(adata)

xy = np.asarray(adata.obsm["X_umap"], dtype=float)
confidence = adata.obs["velocity_confidence"].to_numpy(dtype=float)
latent_time = adata.obs["latent_time"].to_numpy(dtype=float)
spliced = adata.obs["initial_size_spliced"].to_numpy(dtype=float)
unspliced = adata.obs["initial_size_unspliced"].to_numpy(dtype=float)
log_spliced = np.log10(spliced + 1)
log_unspliced = np.log10(unspliced + 1)

# Save cell and gene metrics
metrics = pd.DataFrame({
    "barcode": adata.obs_names,
    "orig.ident": adata.obs["orig.ident"].astype(str).to_numpy(),
    "cell_type": adata.obs["cell_type"].astype(str).to_numpy(),
    "velocity_confidence": confidence,
    "velocity_confidence_transition": adata.obs["velocity_confidence_transition"].to_numpy(dtype=float),
    "latent_time": latent_time,
    "initial_size_spliced": spliced,
    "initial_size_unspliced": unspliced,
    "unspliced_fraction": unspliced / (spliced + unspliced),
})
metrics.to_csv(TABLES / "cell_metrics.csv", index=False)

gene_fit = adata.var[
    ["fit_likelihood", "fit_r2", "fit_alpha", "fit_beta", "fit_gamma", "velocity_genes"]
].copy()
gene_fit.insert(0, "gene", gene_fit.index.astype(str))
gene_fit.to_csv(TABLES / "gene_fit.csv", index=False)

# Plot individual metrics
individual = [
    (confidence, "Velocity confidence", "confidence", "magma", 0, 1, "confidence_umap"),
    (latent_time, "Latent time", "latent time", "viridis", 0, 1, "latent_time_umap"),
    (log_spliced, "Spliced signal", "log10(spliced + 1)", "Blues", None, None, "spliced_umap"),
    (log_unspliced, "Unspliced signal", "log10(unspliced + 1)", "Reds", None, None, "unspliced_umap"),
]
for values, title, label, cmap, vmin, vmax, stem in individual:
    fig, ax = plt.subplots(figsize=(85 / 25.4, 75 / 25.4))
    umap_panel(ax, xy, values, title, label, cmap, vmin, vmax)
    fig.subplots_adjust(left=0.13, right=0.87, bottom=0.13, top=0.90)
    save_figure(fig, stem)

fit_values = pd.to_numeric(gene_fit["fit_likelihood"], errors="coerce").dropna().to_numpy()
fig, ax = plt.subplots(figsize=(85 / 25.4, 75 / 25.4))
ax.hist(fit_values, bins=35, color="#4C72B0", edgecolor="white", linewidth=0.35)
ax.axvline(np.median(fit_values), color="#C44E52", linewidth=0.8, linestyle="--", label=f"median = {np.median(fit_values):.3f}")
ax.set_title("Dynamical fit likelihood")
ax.set_xlabel("fit likelihood (genes)")
ax.set_ylabel("Number of genes")
ax.legend(frameon=False)
ax.tick_params(length=2, width=0.4)
for spine in ax.spines.values():
    spine.set_linewidth(0.4)
fig.subplots_adjust(left=0.16, right=0.96, bottom=0.17, top=0.90)
save_figure(fig, "fit_likelihood")

# Plot overview
fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 120 / 25.4))
umap_panel(axes[0, 0], xy, confidence, "Velocity confidence", "confidence", "magma", 0, 1)
umap_panel(axes[0, 1], xy, latent_time, "Latent time", "latent time", "viridis", 0, 1)
umap_panel(axes[1, 0], xy, log_spliced, "Spliced signal", "log10(spliced + 1)", "Blues")
umap_panel(axes[1, 1], xy, log_unspliced, "Unspliced signal", "log10(unspliced + 1)", "Reds")
axes[0, 2].hist(fit_values, bins=35, color="#4C72B0", edgecolor="white", linewidth=0.3)
axes[0, 2].axvline(np.median(fit_values), color="#C44E52", linewidth=0.8, linestyle="--")
axes[0, 2].set_title("Dynamical fit likelihood", pad=3)
axes[0, 2].set_xlabel("fit likelihood (genes)")
axes[0, 2].set_ylabel("Number of genes")
axes[0, 2].tick_params(length=2, width=0.4)
axes[1, 2].axis("off")
axes[1, 2].text(
    0, 0.95,
    f"Cells: {adata.n_obs:,}\nFitted genes: {len(fit_values):,}\n"
    f"Median confidence: {np.nanmedian(confidence):.3f}\n"
    f"Median latent time: {np.nanmedian(latent_time):.3f}\n"
    f"Median fit likelihood: {np.median(fit_values):.3f}",
    ha="left", va="top", fontsize=8, linespacing=1.5,
)
fig.subplots_adjust(left=0.07, right=0.97, bottom=0.09, top=0.95, wspace=0.38, hspace=0.34)
save_figure(fig, "overview")
