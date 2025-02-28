source("utils/libraries.R")
source("utils/filepaths.R")
library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

# load data
sce <- readRDS(paste0("path_to_GSM8528725_MPS.RDS"))
sce <- scuttle::addPerCellQC(sce)
crop_multiome <- readRDS("path_to_CROP-Multiome.RDS")

# import data from other technologies
qc_data <- read.csv(paste0("path_to_qc_Source_Data.csv"), header = TRUE)

# pull data in same format from sce
add_sce_qc_data <- function(sce, qc_data, modality, nCounts_col, nFeatures_col, method = NA, cell_line = NA) {
  # Extract QC data based on modality
  sce_qc <- sce |> 
    colData() |> 
    as.data.frame() |> 
    rownames_to_column("cell") |>
    dplyr::select(
      nCounts = all_of(nCounts_col), 
      nFeatures = all_of(nFeatures_col), 
      cell = "cell"
    ) |>
    mutate(
      Method = method, 
      Cell_line = cell_line,
      Modality = modality
    )
  
  # Combine and return updated QC data
  qc_data <- rbind(qc_data, sce_qc)
  return(qc_data)
}

# ---------------------------------------------------------------------
# Add datasets from different modalities to qc_data
qc_data <- add_sce_qc_data(sce, qc_data, "RNA", "sum", "detected", method = "MultiPerturb-seq")
qc_data <- add_sce_qc_data(sce, qc_data, "ATAC", "altexps_ATAC_sum", "altexps_ATAC_detected", method = "MultiPerturb-seq")

# Add CROP-seq multiome data
crop_qc_rna <- data.frame(
  nCounts = Matrix::colSums(counts(crop_multiome)),
  nFeatures = Matrix::colSums(counts(crop_multiome) > 0),
  cell = colnames(crop_multiome),
  Method = "CROP-seq multiome",
  Cell_line = NA,
  Modality = "RNA"
)

crop_qc_atac <- data.frame(
  nCounts = Matrix::colSums(counts(altExp(crop_multiome, "ATAC"))),
  nFeatures = Matrix::colSums(counts(altExp(crop_multiome, "ATAC")) > 0),
  cell = colnames(crop_multiome),
  Method = "CROP-seq multiome",
  Cell_line = NA,
  Modality = "ATAC"
)

# Add Seurat object data
seurat_obj <- readRDS("path_to_07_end_integrated_sobj.rds")

seurat_qc <- data.frame(
  nCounts = seurat_obj$nCount_RNA,
  nFeatures = seurat_obj$nFeature_RNA,
  cell = colnames(seurat_obj),
  Method = "CAT-ATAC",
  Cell_line = NA,
  Modality = "RNA"
)

seurat_atac_qc <- data.frame(
  nCounts = seurat_obj$nCount_ATAC,
  nFeatures = seurat_obj$nFeature_ATAC,
  cell = colnames(seurat_obj),
  Method = "CAT-ATAC",
  Cell_line = NA,
  Modality = "ATAC"
)

# Combine all QC data
qc_data <- rbind(qc_data, seurat_qc, seurat_atac_qc, crop_qc_rna, crop_qc_atac)

# Define method order
method_order <- c(
  base::setdiff(unique(qc_data$Method), c("CROP-seq multiome", "MultiPerturb-seq", "CAT-ATAC")),
  "CROP-seq multiome", "MultiPerturb-seq", "CAT-ATAC"
)

# ---------------------------------------------------------------------
# Define plotting parameters and themes
# ---------------------------------------------------------------------

# Define a consistent color palette for methods
method_colors <- c(
  "CROP-seq multiome" = "#E64B35",
  "MultiPerturb-seq" = "#4DBBD5", 
  "CAT-ATAC" = "#00A087",
  "SHARE-seq" = "#3C5488",
  "SNARE-seq" = "#F39B7F",
  "Paired-seq" = "#FF9F00",
  "sci-CAR-seq" = "#8491B4"
)

# Define modality colors
modality_colors <- c(
  "ATAC" = "#3288bd",
  "RNA" = "#d53e4f",
  "CRISPR" = "#fee08b"
)

# Define consistent theme
plot_theme <- theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(size = 14, face = "bold"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    legend.position = "none",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9)
  )

# ---------------------------------------------------------------------
# Filter and prepare data for plotting
# ---------------------------------------------------------------------

# Filter RNA data
qc_data_rna <- qc_data |> 
  filter(Modality == "RNA") |> 
  filter(nCounts > 1000)

qc_data_rna$Method <- factor(qc_data_rna$Method, levels = method_order)

# Prepare capture types for RNA
capture_types_rna <- qc_data |> 
  group_by(Method) |> 
  select(Method, Modality) |> 
  unique() |> 
  filter(!is.na(Method)) |> 
  mutate(
    Modality = factor(Modality, levels = c("CRISPR", "RNA", "ATAC"))
  ) |> 
  rbind(data.frame(Method = "MultiPerturb-seq", Modality = "CRISPR")) |>
  rbind(data.frame(Method = "CAT-ATAC", Modality = "CRISPR")) |>
  rbind(data.frame(Method = "CROP-seq multiome", Modality = "CRISPR")) |>
  mutate(Method = factor(Method, levels = method_order))

# Filter ATAC data
qc_data_atac <- qc_data |> 
  filter(Modality == "ATAC") |> 
  filter(nCounts > 1000)

qc_data_atac$Method <- factor(qc_data_atac$Method, levels = method_order)

# Prepare capture types for ATAC
capture_types_atac <- qc_data |> 
  group_by(Method) |> 
  select(Method, Modality) |> 
  unique() |> 
  rbind(data.frame(Method = "MultiPerturb-seq", Modality = "CRISPR")) |>
  rbind(data.frame(Method = "CAT-ATAC", Modality = "CRISPR")) |>
  rbind(data.frame(Method = "CROP-seq multiome", Modality = "CRISPR")) |>
  mutate(
    Method = factor(Method, levels = method_order),
    Modality = factor(Modality, levels = c("CRISPR", "RNA", "ATAC"))
  ) |> 
  filter(!is.na(Method))

# ---------------------------------------------------------------------
# Create plots
# ---------------------------------------------------------------------

# RNA fragments plot
p_rna_frags <- qc_data_rna %>%
  filter(Modality == "RNA") %>%
  ggplot(aes(x = Method, y = nCounts, fill = Method)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0.5, 
              outlier.alpha = 0.3, notch = TRUE, linewidth = 0.3) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(trans = 'log1p', 
                    breaks = c(1, 10, 100, 1000, 10000, 100000),
                    labels = scales::comma) +
  labs(
    x = "",
    y = "RNA fragments per cell",
    title = "RNA Fragment Distribution"
  ) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(r = 10)))

# ATAC fragments plot
p_atac_frags <- qc_data_atac %>%
  filter(Modality == "ATAC") %>%
  ggplot(aes(x = Method, y = nCounts, fill = Method)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0.5,
              outlier.alpha = 0.3, notch = TRUE, linewidth = 0.3) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(trans = 'log1p',
                    breaks = c(1, 10, 100, 1000, 10000, 50000, 200000),
                    labels = scales::comma) +
  labs(
    x = "Method",
    y = "ATAC fragments per cell",
    title = "ATAC Fragment Distribution"
  ) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(r = 10)))

qc_data_rna <- qc_data_rna |> 
  filter(Modality == "RNA") |> 
  filter(nFeatures > 1000)

qc_data_atac <- qc_data_atac |> 
  filter(Modality == "ATAC") |> 
  filter(nFeatures > 1000)

# RNA features plot
p_rna_feat <- qc_data_rna %>%
  filter(Modality == "RNA") %>%
  ggplot(aes(x = Method, y = nFeatures, fill = Method)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0.5,
              outlier.alpha = 0.3, notch = FALSE, linewidth = 0.3) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(trans = 'log1p',
                    breaks = c(1000, 5000, 10000, 50000),
                    labels = scales::comma) +
  labs(
    x = "",
    y = "Genes per cell",
    title = "Gene Distribution"
  ) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(r = 10)))

# ATAC features plot
p_atac_feat <- qc_data_atac %>%
  filter(Modality == "ATAC") %>%
  ggplot(aes(x = Method, y = nFeatures, fill = Method)) +
  geom_violin(alpha = 0.7, linewidth = 0.3, scale = "width") +
  geom_boxplot(width = 0.2, fill = "white", color = "black", outlier.size = 0.5,
              outlier.alpha = 0.3, notch = FALSE, linewidth = 0.3) +
  scale_fill_manual(values = method_colors) +
  scale_y_continuous(trans = 'log1p',
                    breaks = c(1, 10, 100, 1000, 10000, 50000, 200000),
                    labels = scales::comma) +
  labs(
    x = "",
    y = "ATAC peaks per cell",
    title = "ATAC Peak Distribution"
  ) +
  plot_theme +
  theme(axis.title.y = element_text(margin = margin(r = 10)))

# Create final combined plots (simplified)
combined_frags <- cowplot::plot_grid(
  p_rna_frags + theme(axis.text.x = element_blank(), axis.title.x = element_blank()), 
  p_atac_frags,
  ncol = 1, 
  nrow = 2,
  align = "v",
  labels = c("A", "B"),
  label_size = 16,
  rel_heights = c(1, 1.2)
)

combined_features <- cowplot::plot_grid(
  p_rna_feat + theme(axis.text.x = element_blank(), axis.title.x = element_blank()),
  p_atac_feat,
  ncol = 1, nrow = 2,
  align = "v",
  labels = c("C", "D"),
  label_size = 16,
  rel_heights = c(0.95, 1)
)

# # Save plots with high resolution
# ggsave("qc_fragments_distribution.pdf", combined_frags, width = 10, height = 12, dpi = 300)
# ggsave("qc_features_distribution.pdf", combined_features, width = 10, height = 12, dpi = 300)

# Display plots in R
combined_frags
combined_features