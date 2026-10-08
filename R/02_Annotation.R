# Setup and Configuration
library(tidyverse)
library(Seurat)
library(BiocParallel)
library(qs2)
library(harmony)
library(patchwork)

output_dir <- file.path(getwd(), "output", "standard_pipeline")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

celltype_cols <- c(
    RG = "#6A51A3",
    GIPC = "#D7301F",
    OPC = "#F5B70A",
    ASC = "#FCAE6B",
    NIPC = "#315BA8",
    EN = "#7FB3D5",
    IN = "#00B8A9",
    CR = "#C44E52",
    MG = "#009E73",
    Endo = "#7F9A3A",
    Epend = "#B5B35C",
    Fibro = "#8E6C8A",
    RBC = "#8C3B3B",
    lowQC = "#9E9E9E"
)

# Setup Plot Function
plot_dimplot <- function(
  obj,
  group,
  reduction = "umap",
  colors = NULL,
  pt_size = 0.2,
  pt_alpha = 1,
  axis_text_pt = 6,
  show_labels = TRUE,
  label_text_pt = 5,
  legend_text_pt = 6,
  legend_key_size_cm = 0.22,
  legend_point_size_mm = 1.0,
  axis_linewidth = 0.35,
  seed = 123,
  raster_dpi = 600
) {
    # 1. Extraction
    coord_names <- colnames(obj[[reduction]])[1:2]
    vars <- unique(c(coord_names, group))
    df <- Seurat::FetchData(obj, vars = vars)
    axis_prefix <- toupper(reduction)

    dim_data <- data.frame(
        DIM_1 = df[[coord_names[1]]],
        DIM_2 = df[[coord_names[2]]],
        group_for_points = df[[group]],
        check.names = FALSE
    )

    if (is.factor(dim_data$group_for_points)) {
        dim_data$group_for_points <- droplevels(dim_data$group_for_points)
    }

    # 2. Axis
    min_x <- min(dim_data$DIM_1, na.rm = TRUE)
    min_y <- min(dim_data$DIM_2, na.rm = TRUE)
    max_x <- max(dim_data$DIM_1, na.rm = TRUE)

    arrow_len <- (max_x - min_x) * 0.15
    offset <- (max_x - min_x) * 0.04
    origin_x <- min_x - offset
    origin_y <- min_y - offset
    label_gap <- (max_x - min_x) * 0.02

    x_axis_label <- paste0(axis_prefix, "_1")
    y_axis_label <- paste0(axis_prefix, "_2")

    # 3. Labels
    if (isTRUE(show_labels)) {
        label_data <- dim_data %>%
            dplyr::group_by(group_for_points) %>%
            dplyr::summarise(
                DIM_1 = median(DIM_1, na.rm = TRUE),
                DIM_2 = median(DIM_2, na.rm = TRUE),
                .groups = "drop"
            )
    }

    # 4. Plot
    p <- ggplot2::ggplot(dim_data, ggplot2::aes(x = DIM_1, y = DIM_2)) +
        ggrastr::rasterise(
            ggplot2::geom_point(
                ggplot2::aes(color = group_for_points),
                size = pt_size,
                alpha = pt_alpha,
                shape = 16,
                stroke = 0
            ),
            dpi = raster_dpi
        )

    if (isTRUE(show_labels)) {
        p <- p +
            ggrepel::geom_label_repel(
                data = label_data,
                ggplot2::aes(
                    x = DIM_1,
                    y = DIM_2,
                    label = group_for_points,
                    color = group_for_points
                ),
                fill = "white",
                label.size = 0.15,
                label.r = grid::unit(0.05, "lines"),
                label.padding = grid::unit(0.1, "lines"),
                size = label_text_pt / ggplot2::.pt,
                fontface = "bold",
                seed = seed,
                max.overlaps = Inf,
                show.legend = FALSE
            )
    }

    # 5. Arrows
    p <- p +
        ggplot2::annotate(
            "segment",
            x = origin_x,
            xend = origin_x + arrow_len,
            y = origin_y,
            yend = origin_y,
            arrow = grid::arrow(type = "open", length = grid::unit(0.045, "inches")),
            linewidth = axis_linewidth,
            color = "black"
        ) +
        ggplot2::annotate(
            "text",
            x = origin_x + arrow_len * 0.05,
            y = origin_y - label_gap,
            label = x_axis_label,
            hjust = 0,
            vjust = 1,
            size = axis_text_pt / ggplot2::.pt,
            fontface = "bold"
        ) +
        ggplot2::annotate(
            "segment",
            x = origin_x,
            xend = origin_x,
            y = origin_y,
            yend = origin_y + arrow_len,
            arrow = grid::arrow(type = "open", length = grid::unit(0.045, "inches")),
            linewidth = axis_linewidth,
            color = "black"
        ) +
        ggplot2::annotate(
            "text",
            x = origin_x - label_gap,
            y = origin_y + arrow_len * 0.05,
            label = y_axis_label,
            hjust = 0,
            vjust = 0,
            angle = 90,
            size = axis_text_pt / ggplot2::.pt,
            fontface = "bold"
        ) +
        ggplot2::theme_void() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(3, 3, 3, 3, unit = "mm"),
            legend.position = "right",
            legend.title = ggplot2::element_blank(),
            legend.text = ggplot2::element_text(size = legend_text_pt, lineheight = 1.5),
            legend.key.size = grid::unit(legend_key_size_cm, "cm"),
            legend.key.height = grid::unit(legend_key_size_cm * 1.5, "cm"),
            legend.spacing.y = grid::unit(0.04, "cm")
        ) +
        ggplot2::coord_fixed(clip = "off")

    # 6. Colors
    if (is.null(colors)) {
        p <- p + ggplot2::scale_color_discrete()
    } else {
        p <- p + ggplot2::scale_color_manual(values = colors)
    }

    p <- p +
        ggplot2::guides(
            color = ggplot2::guide_legend(
                ncol = 1,
                byrow = TRUE,
                override.aes = list(
                    size = legend_point_size_mm,
                    alpha = 1,
                    shape = 16,
                    stroke = 0
                )
            )
        )

    return(p)
}

# Load data
TNC <- qs_read("TNC.qs2")

# Cell type annotation
plot_dimplot(
    TNC,
    group = "seurat_clusters",
    reduction = "umap"
)

Idents(TNC) <- TNC$seurat_clusters

TNC$cell_type <- "Undefined"

# Radial glia
FeaturePlot(TNC,
    features = c(
        "Pax6", "Sox2", "Fabp7",
        "Hes1", "Hes5", "Vim",
        "Nes", "Slc1a3", "Ptprz1"
    ),
    order = TRUE,
    ncol = 3,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(2)
    )
] <- "RG"

# Neuronal intermediate progenitor cells
FeaturePlot(TNC,
    features = c(
        "Eomes", "Neurog2", "Neurod1",
        "Neurod2", "Neurod6"
    ),
    order = TRUE,
    ncol = 3,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(11)
    )
] <- "NIPC"

# Glial intermediate progenitor cells
FeaturePlot(TNC,
    features = c(
        "Ascl1", "Egfr",
        "Sox9", "Dbi"
    ),
    order = TRUE,
    ncol = 2,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(6)
    )
] <- "GIPC"

# Oligodendrocyte precursor cells
FeaturePlot(TNC,
    features = c("Pdgfra", "Sox10"),
    order = TRUE,
    ncol = 2,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(9)
    )
] <- "OPC"

# Astrocytes
FeaturePlot(TNC,
    features = c(
        "Aldh1l1", "Aqp4", "Gfap",
        "Slc1a2", "Slc1a3", "Aldoc",
        "S100b", "Gja1", "Sparcl1"
    ),
    order = TRUE,
    ncol = 3,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(4)
    )
] <- "ASC"

# Excitatory Neurons
FeaturePlot(TNC,
    features = c(
        "Slc17a7", "Slc17a6", "Neurod2",
        "Neurod6", "Tbr1", "Satb2",
        "Bcl11b", "Cux1", "Cux2"
    ),
    order = TRUE,
    ncol = 3,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(0, 1, 5, 13, 17)
    )
] <- "EN"

# Inhibitory Neurons
FeaturePlot(TNC,
    features = c(
        "Gad1", "Gad2", "Slc32a1",
        "Dlx1", "Dlx2", "Dlx5",
        "Dlx6", "Lhx6", "Adarb2"
    ),
    order = TRUE,
    ncol = 3,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(3, 8, 10, 12, 14, 15, 16)
    )
] <- "IN"

# Cajal-Retzius cells
FeaturePlot(TNC,
    features = c("Reln"),
    order = TRUE,
    ncol = 1,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(20)
    )
] <- "CR"

# Microglia
FeaturePlot(TNC,
    features = c("Cx3cr1", "Aif1"),
    order = TRUE,
    ncol = 2,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(19)
    )
] <- "MG"

# Endothelial cells
FeaturePlot(TNC,
    features = c("Cldn5", "Pecam1"),
    order = TRUE,
    ncol = 2,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(18)
    )
] <- "Endo"

# Fibroblasts
FeaturePlot(TNC,
    features = c("Col1a1", "Dcn", "Lum"),
    order = TRUE,
    ncol = 3,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(21)
    )
] <- "Fibro"

# Ependymal cells
FeaturePlot(TNC,
    features = c("Foxj1", "Ccdc153"),
    order = TRUE,
    ncol = 2,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(22)
    )
] <- "Epend"

# Red blood cells
FeaturePlot(TNC,
    features = c("Hba-a1"),
    order = TRUE,
    ncol = 1,
    cols = c("lightgrey", "red")
) & coord_fixed()

TNC$cell_type[
    WhichCells(TNC,
        ident = c(23)
    )
] <- "RBC"

# Low-quality cells
TNC[["percent.ribo"]] <- PercentageFeatureSet(
    TNC,
    pattern = "^Rp[sl]"
)

FeaturePlot(
    TNC,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    reduction = "umap",
    ncol = 2
)

TNC$cell_type[
    WhichCells(TNC,
        ident = c(7)
    )
] <- "lowQC"

# Annotation visualization
p <- plot_dimplot(
    TNC,
    group = "cell_type",
    reduction = "umap",
    colors = celltype_cols
)

ggsave(
    filename = file.path(output_dir, "UMAP_total.pdf"),
    plot = p,
    width = 9,
    height = 9,
    units = "cm"
)

# Save Point
qs_save(TNC, "TNC.qs2")
