# visium_smoking_dataset_seurat_workflow.R
# Purpose:
#   - Load multiple 10x Visium lung samples (Current/Past/Never smokers)
#   - Perform basic normalisation + QC + merging
#   - Run PCA/UMAP/clustering
#   - Score a smoking gene set with AddModuleScore and visualise on spatial images
#
# Data provenance:
#   - These Visium datasets were collected from public datasets provided in:
#       PMID: 36543915
#
# Repo assumptions (relative to repo root):
#   - data/visium/WSA_LngSP10193347_A37_CS/
#   - data/visium/WSA_LngSP9258464_A37_CS/
#   - data/visium/WSA_LngSP9258469_A47_NS/
#   - data/visium/WSA_LngSP9258465_A47_NS/
#   - data/visium/WSA_LngSP9258463_A50_PS/
#   - data/visium/WSA_LngSP9258467_A50_PS/
#   - A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke/Alen_smoke_geneset.csv
#
# Outputs:
#   - Figures saved to results/figures/
#

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
})

# -----------------------------
# User-configurable parameters
# -----------------------------
# Visium sample directories (edit to match your repo structure)
visium_root <- file.path("data", "visium")

sample_map <- list(
  CS1 = list(dir = file.path(visium_root, "WSA_LngSP10193347_A37_CS"), slice = "Current.smoker.1", group = "CS", mito_max = 10),
  CS2 = list(dir = file.path(visium_root, "WSA_LngSP9258464_A37_CS"),  slice = "Current.smoker.2", group = "CS", mito_max = 10),
  NS1 = list(dir = file.path(visium_root, "WSA_LngSP9258469_A47_NS"),  slice = "Never.smoker.1",   group = "NS", mito_max = 20),
  NS2 = list(dir = file.path(visium_root, "WSA_LngSP9258465_A47_NS"),  slice = "Never.smoker.2",   group = "NS", mito_max = 20),
  PS1 = list(dir = file.path(visium_root, "WSA_LngSP9258463_A50_PS"),  slice = "Past.smoker.1",    group = "PS", mito_max = 10),
  PS2 = list(dir = file.path(visium_root, "WSA_LngSP9258467_A50_PS"),  slice = "Past.smoker.2",    group = "PS", mito_max = 10)
)

# Gene set file
geneset_dir <- "A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke"
geneset_file <- file.path(geneset_dir, "Alen_smoke_geneset.csv")

# Plot/output settings
plot_dir <- file.path("results", "figures")
save_plots <- TRUE

# Dimensionality reduction parameters
n_pcs <- 30
cluster_resolution <- 0.1

# Module score visualisation (continuous scale)
module_score_limits <- c(0, 0.2)
module_score_breaks <- c(0, 0.05, 0.1, 0.15, 0.2)

# -----------------------------
# Checks + setup
# -----------------------------
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(geneset_dir)) {
  stop(
    "Directory not found: ", geneset_dir, "\n",
    "Tip: set your working directory to the GitHub repository root before running this script."
  )
}
if (!file.exists(geneset_file)) stop("Missing gene set file: ", geneset_file)

# Optional dependency note (don’t auto-install inside scripts)
if (utils::packageVersion("Matrix") < "1.5.3") {
  message(
    "NOTE: Matrix version is ", utils::packageVersion("Matrix"),
    ". If you hit spatial dependency issues, consider updating Matrix."
  )
}

# -----------------------------
# Helper functions
# -----------------------------
read_visium_sample <- function(data_dir, slice_name, assay_name = "Spatial", image_name = "tissue_lowres_image.png") {
  spatial_dir <- file.path(data_dir, "spatial")
  
  if (!dir.exists(data_dir)) stop("Sample directory not found: ", data_dir)
  if (!dir.exists(spatial_dir)) stop("Missing 'spatial' folder in: ", data_dir)
  
  img <- Read10X_Image(
    image.dir = spatial_dir,
    image.name = image_name
  )
  
  obj <- Load10X_Spatial(
    data.dir = data_dir,
    filename = "filtered_feature_bc_matrix.h5",
    filter.matrix = TRUE,
    to.upper = FALSE,
    assay = assay_name,
    slice = slice_name,
    image = img
  )
  obj
}

save_gg <- function(p, filename, width = 6, height = 6, dpi = 600, bg = "white") {
  if (!isTRUE(save_plots)) return(invisible(NULL))
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = bg
  )
}

# -----------------------------
# Load samples
# -----------------------------
message("Loading Visium samples (PMID: 36543915 public datasets)...")

objs <- list()
for (nm in names(sample_map)) {
  info <- sample_map[[nm]]
  message("  - ", nm, " from: ", info$dir)
  
  obj <- read_visium_sample(
    data_dir = info$dir,
    slice_name = info$slice
  )
  
  # Add sample metadata
  obj$sample <- nm
  obj$smoking_group <- info$group
  
  objs[[nm]] <- obj
}

# Quick QC plots: nCount_Spatial on each sample
for (nm in names(objs)) {
  p <- SpatialFeaturePlot(
    objs[[nm]],
    features = "nCount_Spatial",
    image.alpha = 1,
    pt.size.factor = 3
  ) + ggtitle(paste0(nm, " - nCount_Spatial"))
  print(p)
}

# -----------------------------
# Normalisation + scaling + mito QC
# -----------------------------
message("Normalising, scaling, and computing percent mitochondria...")

for (nm in names(objs)) {
  obj <- objs[[nm]]
  
  obj <- NormalizeData(obj, assay = "Spatial", verbose = TRUE)
  obj <- ScaleData(obj, assay = "Spatial", verbose = FALSE)
  
  # Mito percent (assuming human-style MT- prefix)
  obj <- PercentageFeatureSet(obj, pattern = "^MT-", col.name = "percent_mito")
  
  # Basic QC violin (printed; you can save if you want)
  v <- VlnPlot(obj, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito"), ncol = 3) +
    ggtitle(paste0(nm, " QC"))
  print(v)
  
  # Filter based on per-sample threshold
  mito_max <- sample_map[[nm]]$mito_max
  obj <- obj[, obj$percent_mito < mito_max]
  
  objs[[nm]] <- obj
}

# -----------------------------
# Merge all samples
# -----------------------------
message("Merging samples...")
ALL <- merge(objs[[1]], y = objs[-1])

# -----------------------------
# Feature selection + PCA/UMAP/Clustering
# -----------------------------
message("Running variable feature selection, PCA, neighbors, clusters, UMAP...")

DefaultAssay(ALL) <- "Spatial"

ALL <- FindVariableFeatures(ALL, selection.method = "vst", nfeatures = 2000)
ALL <- RunPCA(ALL, verbose = TRUE)

ALL <- FindNeighbors(ALL, reduction = "pca", dims = 1:n_pcs)
ALL <- FindClusters(ALL, resolution = cluster_resolution)

# NOTE: your original script used assay = "SCT" without SCTransform.
# Here we run UMAP on the existing PCA from the Spatial assay.
ALL <- RunUMAP(ALL, reduction = "pca", dims = 1:n_pcs)

# Check integration/merge quality
p_by_sample <- DimPlot(ALL, reduction = "umap", group.by = "sample")
p_by_cluster <- DimPlot(ALL, reduction = "umap", label = FALSE, pt.size = 1)

print(p_by_sample)
print(p_by_cluster)

save_gg(p_by_sample, "UMAP_by_sample.png", width = 7, height = 5, dpi = 400)
save_gg(p_by_cluster, "UMAP_by_cluster.png", width = 7, height = 5, dpi = 400)

# Spatial overview (shows all images)
p_spatial <- SpatialDimPlot(ALL, image.alpha = 1, pt.size.factor = 3)
print(p_spatial)
save_gg(p_spatial, "SpatialDimPlot_all_samples.png", width = 10, height = 7, dpi = 400)

# -----------------------------
# Load smoking gene set and score
# -----------------------------
message("Loading gene set and computing module scores...")

s <- read.csv(geneset_file, stringsAsFactors = FALSE)

required_cols <- c("log2FoldChange", "padj", "hgnc_symbol")
missing_cols <- setdiff(required_cols, colnames(s))
if (length(missing_cols) > 0) {
  stop("Gene set file missing required columns: ", paste(missing_cols, collapse = ", "))
}

geneset <- s %>%
  filter(log2FoldChange > 0, padj < 0.05) %>%
  pull(hgnc_symbol) %>%
  na.omit() %>%
  unique()

if (length(geneset) < 20) {
  message("NOTE: Gene set has only ", length(geneset), " genes after filtering.")
}

top10 <- head(geneset, 10)
top20 <- head(geneset, 20)

present <- intersect(geneset, rownames(ALL[["Spatial"]]))
message("Genes present in ALL (Spatial assay): ", length(present), " / ", length(geneset))

# Seurat v5 compatibility: JoinLayers if available/needed
if ("JoinLayers" %in% getNamespaceExports("Seurat")) {
  ALL <- JoinLayers(ALL)
}

# AddModuleScore creates columns with suffix "1"
ALL <- AddModuleScore(ALL, features = list(geneset), name = "smoking_gene_set", assay = "Spatial")
ALL <- AddModuleScore(ALL, features = list(top10),  name = "smoking_gene_set_top10", assay = "Spatial")
ALL <- AddModuleScore(ALL, features = list(top20),  name = "smoking_gene_set_top20", assay = "Spatial")

# Detect the actual column names created (avoid fragile index-based renaming)
score_cols <- c(
  smoking_gene_set = grep("^smoking_gene_set\\d+$", colnames(ALL@meta.data), value = TRUE),
  smoking_gene_set_top10 = grep("^smoking_gene_set_top10\\d+$", colnames(ALL@meta.data), value = TRUE),
  smoking_gene_set_top20 = grep("^smoking_gene_set_top20\\d+$", colnames(ALL@meta.data), value = TRUE)
)

if (any(vapply(score_cols, length, integer(1)) != 1)) {
  stop("Unexpected module score columns. Found: ", paste(unlist(score_cols), collapse = ", "))
}

# For plotting convenience, alias to clean names (copies values)
ALL$smoking_gene_set <- ALL[[score_cols["smoking_gene_set"]]][, 1]
ALL$smoking_gene_set_top10 <- ALL[[score_cols["smoking_gene_set_top10"]]][, 1]
ALL$smoking_gene_set_top20 <- ALL[[score_cols["smoking_gene_set_top20"]]][, 1]

# -----------------------------
# Plot module score on a specific image
# -----------------------------
# Example: plot on "Never.smoker.1" image (NS1 slice name)
image_to_plot <- "Never.smoker.1"

p_score <- SpatialFeaturePlot(
  ALL,
  features = "smoking_gene_set",
  alpha = 1,
  images = image_to_plot
) +
  scale_fill_continuous(
    limits = module_score_limits,
    breaks = module_score_breaks
  ) +
  ggtitle(paste0("Smoking gene set score - ", image_to_plot, " (PMID: 36543915)"))

print(p_score)
save_gg(p_score, "smoking_gene_set_NS1.png", width = 6, height = 6, dpi = 1000)

message("Done. Outputs saved to: ", plot_dir)
