# hlca_module_scores.R
# Purpose:
#   - Load two gene sets (smoking signature + primate cluster 2)
#   - Filter to genes present in HLCA Seurat object
#   - Score each gene set using Seurat::AddModuleScore
#   - Visualize scores on HLCA UMAP
#
# Assumptions:
#   - HLCA Seurat object is already loaded in the environment as `HLCA`; available from official website --- The integrated Human Lung Cell Atlas (HLCA) v1.0
#   - Your GitHub repo contains this folder at the repo root:
#       A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke/
#         ├─ cluster_2_genes.csv
#         └─ de_genelist_smoking_NORM_annotated.csv
#
# How to run (from repo root):
#   source("scripts/hlca_module_scores.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# User-configurable parameters
# -----------------------------
geneset_dir <- "A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke"

smoking_file <- file.path(geneset_dir, "de_genelist_smoking_NORM_annotated.csv")
primate_file <- file.path(geneset_dir, "cluster_2_genes.csv")

padj_cutoff <- 0.05
log2fc_cutoff <- 1

reduction_name <- "umap"
pt_size <- 0.1

# Optional: save plots
save_plots <- FALSE
plot_dir <- file.path("results", "figures")
plot_width <- 7
plot_height <- 6
plot_dpi <- 300

# -----------------------------
# Basic checks
# -----------------------------
if (!exists("HLCA")) {
  stop("Object `HLCA` not found. Please load the HLCA Seurat object before running this script.")
}

if (!inherits(HLCA, "Seurat")) {
  stop("`HLCA` exists but is not a Seurat object.")
}

# Check expected folder structure (repo root)
if (!dir.exists(geneset_dir)) {
  stop(
    "Directory not found: ", geneset_dir, "\n",
    "Tip: set your working directory to the GitHub repository root before running this script."
  )
}

if (!file.exists(smoking_file)) stop("Missing file: ", smoking_file)
if (!file.exists(primate_file)) stop("Missing file: ", primate_file)

if (!reduction_name %in% Reductions(HLCA)) {
  stop(sprintf(
    "Reduction '%s' not found in HLCA. Available reductions: %s",
    reduction_name, paste(Reductions(HLCA), collapse = ", ")
  ))
}

# -----------------------------
# Load gene sets
# -----------------------------
smk_df <- read.csv(smoking_file, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("padj", "log2FoldChange", "hgnc_symbol")
missing_cols <- setdiff(required_cols, colnames(smk_df))
if (length(missing_cols) > 0) {
  stop(
    "Smoking signature file is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

smk_sig <- smk_df %>%
  filter(padj < padj_cutoff, log2FoldChange > log2fc_cutoff) %>%
  pull(hgnc_symbol) %>%
  na.omit() %>%
  unique()

primate_df <- read.delim(primate_file, header = TRUE, stringsAsFactors = FALSE)

# Your file appears to store genes in a column named `x`
if (!"x" %in% colnames(primate_df)) {
  stop(
    "Primate cluster file does not contain a column named 'x'. Columns found: ",
    paste(colnames(primate_df), collapse = ", ")
  )
}

primate_cluster_2 <- primate_df %>%
  pull(x) %>%
  na.omit() %>%
  unique()

gene_sets <- list(
  smk_sig = smk_sig,
  primate_cluster_2 = primate_cluster_2
)

# -----------------------------
# Helper: keep only genes present
# -----------------------------
filter_present_genes <- function(seu, genes) {
  genes <- unique(stats::na.omit(genes))
  genes[genes %in% rownames(seu)]
}

gene_sets_present <- lapply(gene_sets, function(gs) filter_present_genes(HLCA, gs))

# Fail fast if any gene set becomes empty
empty_sets <- names(gene_sets_present)[vapply(gene_sets_present, length, integer(1)) == 0]
if (length(empty_sets) > 0) {
  stop(
    "These gene sets have 0 genes present in HLCA: ",
    paste(empty_sets, collapse = ", "),
    "\nCheck gene identifiers (HGNC symbols) and HLCA rownames."
  )
}

message("Genes retained after filtering to HLCA features:")
for (nm in names(gene_sets_present)) {
  message(sprintf("  - %s: %d genes", nm, length(gene_sets_present[[nm]])))
}

# -----------------------------
# Score gene sets
# -----------------------------
# AddModuleScore will create metadata columns: <name>1, <name>2, ...
HLCA <- AddModuleScore(
  object = HLCA,
  features = gene_sets_present,
  name = names(gene_sets_present)
)

score_cols <- setNames(paste0(names(gene_sets_present), "1"), names(gene_sets_present))

# -----------------------------
# Plot
# -----------------------------
p_smk <- FeaturePlot(
  object = HLCA,
  features = unname(score_cols["smk_sig"]),
  reduction = reduction_name,
  pt.size = pt_size,
  order = TRUE
) + ggtitle("Module score: smoking signature")

p_pri <- FeaturePlot(
  object = HLCA,
  features = unname(score_cols["primate_cluster_2"]),
  reduction = reduction_name,
  pt.size = pt_size,
  order = TRUE
) + ggtitle("Module score: primate_cluster_2")

print(p_smk)
print(p_pri)

# -----------------------------
# Optional: save plots
# -----------------------------
if (isTRUE(save_plots)) {
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggsave(
    filename = file.path(plot_dir, "HLCA_moduleScore_smoking_signature.png"),
    plot = p_smk,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )
  
  ggsave(
    filename = file.path(plot_dir, "HLCA_moduleScore_primate_cluster_2.png"),
    plot = p_pri,
    width = plot_width,
    height = plot_height,
    dpi = plot_dpi
  )
}

