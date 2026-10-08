# Setup and Configuration
library(tidyverse)
library(Seurat)
library(BiocParallel)
library(qs2)
library(harmony)

output_dir <- file.path(getwd(), "output", "standard_pipeline")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(123)

dataset_dir <- "data"

# Controls
Ctrl1.mat <- Read10X(data.dir = file.path(dataset_dir, "Ctrl1"))
Ctrl1 <- CreateSeuratObject(
	counts = Ctrl1.mat,
	project = "Ctrl1",
	min.cells = 3,
	min.features = 200
)

Ctrl2.mat <- Read10X(data.dir = file.path(dataset_dir, "Ctrl2"))
Ctrl2 <- CreateSeuratObject(
	counts = Ctrl2.mat,
	project = "Ctrl2",
	min.cells = 3,
	min.features = 200
)

Ctrl3.mat <- Read10X(data.dir = file.path(dataset_dir, "Ctrl3"))
Ctrl3 <- CreateSeuratObject(
	counts = Ctrl3.mat,
	project = "Ctrl3",
	min.cells = 3,
	min.features = 200
)

# KOs
KO1.mat <- Read10X(data.dir = file.path(dataset_dir, "KO1"))
KO1 <- CreateSeuratObject(
	counts = KO1.mat,
	project = "KO1",
	min.cells = 3,
	min.features = 200
)

KO2.mat <- Read10X(data.dir = file.path(dataset_dir, "KO2"))
KO2 <- CreateSeuratObject(
	counts = KO2.mat,
	project = "KO2",
	min.cells = 3,
	min.features = 200
)

KO3.mat <- Read10X(data.dir = file.path(dataset_dir, "KO3"))
KO3 <- CreateSeuratObject(
	counts = KO3.mat,
	project = "KO3",
	min.cells = 3,
	min.features = 200
)

rm(list = ls(pattern = "\\.mat$"))

# QC
max_nFeature_RNA <- 7500
max_percent_mt <- 3
plot_qc <- TRUE

qc_report_file <- file.path(output_dir, "TNC_QC_before_after_report.pdf")
dir.create(dirname(qc_report_file), recursive = TRUE, showWarnings = FALSE)

while (!is.null(dev.list())) {
  dev.off()
}

sample_list <- list(
  Ctrl1 = Ctrl1,
  Ctrl2 = Ctrl2,
  Ctrl3 = Ctrl3,
  KO1 = KO1,
  KO2 = KO2,
  KO3 = KO3
)

process_qc <- function(sobj, sample_name) {
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^mt-")
  sobj[["all"]] <- "all"

  if (plot_qc) {
    p1 <- VlnPlot(
      sobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3,
      group.by = "all",
      layer = "counts"
    ) +
      patchwork::plot_annotation(
        title = "Pre-filter"
      )
  }

  meta_df <- sobj[[]]

  keep_cells <- rownames(meta_df)[
    meta_df$nFeature_RNA < max_nFeature_RNA &
      meta_df$percent.mt < max_percent_mt
  ]

  sobj <- subset(sobj, cells = keep_cells)

  if (plot_qc) {
    p2 <- VlnPlot(
      sobj,
      features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
      ncol = 3,
      group.by = "all",
      layer = "counts"
    ) +
      patchwork::plot_annotation(
        title = "Post-filter"
      )

    page_plot <- patchwork::wrap_plots(p1, p2, ncol = 1) +
      patchwork::plot_annotation(
        title = sample_name,
        subtitle = paste0(
          "QC filtering: nFeature_RNA < ", max_nFeature_RNA,
          ", percent.mt < ", max_percent_mt
        ),
        theme = ggplot2::theme(
          plot.title = ggplot2::element_text(
            size = 20,
            face = "bold",
            hjust = 0.5
          ),
          plot.subtitle = ggplot2::element_text(
            size = 12,
            hjust = 0.5
          )
        )
      )

    print(page_plot)
  }

  sobj[["all"]] <- NULL

  return(sobj)
}

grDevices::pdf(qc_report_file, width = 12, height = 10)

tryCatch(
  {
    sample_list <- Map(process_qc, sample_list, names(sample_list))
  },
  finally = {
    grDevices::dev.off()
  }
)

Ctrl1 <- sample_list[["Ctrl1"]]
Ctrl2 <- sample_list[["Ctrl2"]]
Ctrl3 <- sample_list[["Ctrl3"]]
KO1 <- sample_list[["KO1"]]
KO2 <- sample_list[["KO2"]]
KO3 <- sample_list[["KO3"]]

rm(sample_list, process_qc, max_nFeature_RNA, max_percent_mt, plot_qc)

# Merge Seurat objects
TNC <- merge(
  x = Ctrl1, 
  y = c(Ctrl2, Ctrl3, KO1, KO2, KO3), 
  add.cell.ids = c("Ctrl1", "Ctrl2", "Ctrl3", "KO1", "KO2", "KO3"),
)

rm(Ctrl1, Ctrl2, Ctrl3, KO1, KO2, KO3)

# Normalization and Cell Cycle Regression
options(future.globals.maxSize = 16 * 1024^3)

TNC <- JoinLayers(TNC, assays = "RNA")
DefaultAssay(TNC) <- "RNA"
TNC <- NormalizeData(
  TNC,
  assay = "RNA",
  normalization.method = "LogNormalize",
  verbose = FALSE
)

s.genes <- stringr::str_to_title(cc.genes.updated.2019$s.genes)
g2m.genes <- stringr::str_to_title(cc.genes.updated.2019$g2m.genes)

TNC <- CellCycleScoring(
  TNC,
  s.features = s.genes,
  g2m.features = g2m.genes,
  assay = "RNA",
  set.ident = TRUE
)

library(glmGamPoi)
TNC <- SCTransform(
  TNC,
  assay = "RNA",
  new.assay.name = "SCT",
  vars.to.regress = c("S.Score", "G2M.Score"),
  variable.features.n = 2000,
  verbose = FALSE
)

# Save Point
qs_save(TNC, "TNC.qs2")
# TNC <- qs_read("TNC.qs2")

# PCA
TNC <- RunPCA(TNC, verbose = FALSE)

# Elbow plot
ElbowPlot(TNC, ndims = 50, reduction = "pca")

# Save PCA loadings
pca_loadings_pdf <- file.path(output_dir, "PCA_loadings.pdf")
dir.create(dirname(pca_loadings_pdf), recursive = TRUE, showWarnings = FALSE)

max_pcs <- ncol(Embeddings(TNC, reduction = "pca"))
dims <- seq_len(min(30, max_pcs))

pdf_dev <- grDevices::dev.cur()

tryCatch(
  {
    grDevices::pdf(pca_loadings_pdf, width = 9, height = 9, onefile = TRUE)

    for (d in dims) {
      p <- VizDimLoadings(
        TNC,
        dims = d,
        reduction = "pca",
        ncol = 1
      )

      print(p)
    }
  },
  error = function(e) {
    message("Failed to write PCA loadings PDF: ", conditionMessage(e))
  },
  finally = {
    while (grDevices::dev.cur() != pdf_dev) {
      grDevices::dev.off()
    }
  }
)

# Set Dims
dims_use <- c(1, 3:4, 6:7, 9, 11:14, 16, 18:26, 28, 30)

# Harmony
TNC <- RunHarmony(TNC,
  reduction = "pca",
  group.by.vars = "orig.ident",
  reduction.save = "harmony",
  dims = dims_use
)

harmony_dims <- seq_along(dims_use)

TNC <- FindNeighbors(TNC,
  reduction = "harmony",
  dims = harmony_dims
)

TNC <- RunUMAP(
  TNC,
  reduction = "harmony",
  dims = harmony_dims,
  n.neighbors = 40,
  min.dist = 0.3,
  spread = 1.2,
  reduction.name = "umap"
)

TNC <- FindClusters(TNC,
  resolution = 0.5
)

# orig.ident visualization
ident_colors <- c(
  "KO1" = "#3f51b5",
  "KO2" = "#5c6bc0",
  "KO3" = "#7986cb",
  "Ctrl1" = "#fee090",
  "Ctrl2" = "#fc8d59",
  "Ctrl3" = "#d73027"
)

# UMAP Plot
umap_ident <- DimPlot(TNC,
  reduction = "umap", group.by = "orig.ident",
  cols = ident_colors,
  pt.size = 0.05,
  shuffle = TRUE
) +
  coord_fixed(ratio = 1) +
  annotate(
    "text",
    x = -Inf, y = -Inf,
    label = paste0("n = ", ncol(TNC)),
    size = 4,
    color = "black",
    fontface = "bold",
    hjust = -0.2,
    vjust = -1.5
  )

print(umap_ident)

ggsave(
  filename = file.path(output_dir, "UMAP_orig_ident.pdf"),
  plot = umap_ident,
  width = 10,
  height = 8,
  device = cairo_pdf
)

# Save Point
qs_save(TNC, "TNC.qs2")
saveRDS(TNC, "TNC_after_integration.rds")
