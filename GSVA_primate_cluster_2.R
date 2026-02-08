# GSVA_primate_cluster_2.R
# Purpose:
#   - Download GEO datasets (GSE30660, GSE82137)
#   - Map microarray probe IDs to HGNC symbols (hgu133plus2)
#   - Load NORM bronchial biopsy bulk RNA-seq (CPM-normalised)
#   - Convert Ensembl IDs to HGNC symbols (biomaRt)
#   - Run GSVA for primate_cluster_2 gene set across all datasets
#   - Generate and save publication-ready boxplots
#
# Repo assumptions (relative to repo root):
#   - A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke/cluster_2_genes.csv
#   - data/NORM_GEO/cpm_normalised.txt               (you can rename paths below if needed)
#   - metadata/NORM_sample_conditions.csv            (optional; see sample_conditions section)
#

suppressPackageStartupMessages({
  library(GSVA)
  library(GEOquery)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(dplyr)
  library(ggplot2)
  library(ggpubr)
  library(biomaRt)
})

# -----------------------------
# User-configurable parameters
# -----------------------------
geneset_dir <- "A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke"
primate_file <- file.path(geneset_dir, "cluster_2_genes.csv")

# NORM input (change this to match your repo structure)
norm_expr_file <- file.path("data", "NORM_GEO", "cpm_normalised.txt")

# Optional: NORM sample conditions (recommended to make script reproducible)
# CSV with columns: sample, Condition  (Condition values: "Smoker" / "Nonsmoker")
norm_conditions_file <- file.path("metadata", "NORM_sample_conditions.csv")

# Output figures
plot_dir <- file.path("results", "figures")
save_plots <- TRUE
plot_dpi <- 600
plot_bg <- "white"

# Dataset-specific knobs
gse82137_drop_cols <- 4:15   # your original code removed columns 4:15
set.seed(1)

# -----------------------------
# Basic checks
# -----------------------------
if (!dir.exists(geneset_dir)) {
  stop(
    "Directory not found: ", geneset_dir, "\n",
    "Tip: set your working directory to the GitHub repository root before running this script."
  )
}
if (!file.exists(primate_file)) stop("Missing file: ", primate_file)

if (!file.exists(norm_expr_file)) {
  stop(
    "Missing NORM expression file: ", norm_expr_file, "\n",
    "Put it in the repo (recommended) or update `norm_expr_file` to your local path."
  )
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper functions
# -----------------------------
probe_to_symbol_matrix <- function(expr_mat) {
  # Map probe IDs to HGNC symbols using hgu133plus2.db, then collapse duplicates by mean.
  stopifnot(is.matrix(expr_mat) || is.data.frame(expr_mat))
  expr_mat <- as.matrix(expr_mat)
  
  symbols <- mapIds(
    hgu133plus2.db,
    keys = rownames(expr_mat),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  keep <- !is.na(symbols)
  expr_mat <- expr_mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  rownames(expr_mat) <- symbols
  
  # Collapse duplicate symbols by mean
  df <- as.data.frame(expr_mat)
  df$Gene <- rownames(df)
  df <- df %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  out <- as.matrix(df[, -1, drop = FALSE])
  rownames(out) <- df$Gene
  out
}

ensembl_to_hgnc_matrix <- function(expr_df) {
  # Input: data.frame with rownames = Ensembl gene IDs, columns = samples (numeric)
  # Output: matrix with rownames = HGNC symbols, collapsed duplicates by mean
  stopifnot(is.data.frame(expr_df) || is.matrix(expr_df))
  expr_df <- as.data.frame(expr_df)
  
  ensembl_ids <- rownames(expr_df)
  if (is.null(ensembl_ids)) stop("NORM expression input must have rownames as Ensembl IDs.")
  
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  conv <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters = "ensembl_gene_id",
    values = ensembl_ids,
    mart = mart
  ) %>%
    filter(!is.na(hgnc_symbol), hgnc_symbol != "")
  
  if (nrow(conv) == 0) stop("biomaRt returned no HGNC symbols. Check Ensembl IDs / internet access.")
  
  expr_df$ensembl_gene_id <- rownames(expr_df)
  merged <- merge(expr_df, conv, by = "ensembl_gene_id")
  
  # Drop ensembl ID column, collapse duplicate HGNC by mean
  merged <- merged %>%
    select(-ensembl_gene_id) %>%
    group_by(hgnc_symbol) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  out <- as.matrix(merged[, -1, drop = FALSE])
  rownames(out) <- merged$hgnc_symbol
  out
}

run_gsva_single_set <- function(expr_mat, gene_set, kcdf = NULL, verbose = TRUE) {
  # gene_set: character vector
  gene_set <- unique(stats::na.omit(gene_set))
  if (length(gene_set) == 0) stop("Gene set is empty.")
  
  gs <- list(primate_cluster_2 = gene_set)
  
  # Filter gene set to those present
  gs[[1]] <- gs[[1]][gs[[1]] %in% rownames(expr_mat)]
  if (length(gs[[1]]) == 0) stop("No genes from gene set are present in the expression matrix.")
  
  if (is.null(kcdf)) {
    GSVA::gsva(expr_mat, gs, method = "gsva", verbose = verbose)
  } else {
    GSVA::gsva(expr_mat, gs, method = "gsva", kcdf = kcdf, verbose = verbose)
  }
}

make_gsva_boxplot <- function(gsva_vec, condition, ylab = "GSVA score", title = NULL) {
  df <- data.frame(
    GSVA_Score = as.numeric(gsva_vec),
    Condition = as.factor(condition)
  )
  
  ggplot(df, aes(x = Condition, y = GSVA_Score, fill = Condition)) +
    geom_boxplot(color = "black", width = 0.6) +
    geom_jitter(width = 0.2, size = 2, color = "black") +
    labs(title = title, x = NULL, y = ylab, fill = NULL) +
    theme_classic()
}

save_plot <- function(p, filename, width = 4, height = 3) {
  ggsave(
    filename = file.path(plot_dir, filename),
    plot = p,
    width = width,
    height = height,
    units = "in",
    dpi = plot_dpi,
    bg = plot_bg
  )
}

# -----------------------------
# Load gene set
# -----------------------------
primate_df <- read.delim(primate_file, header = TRUE, stringsAsFactors = FALSE)

# Your original file uses a column called `x`
if (!"x" %in% colnames(primate_df)) {
  stop("Gene set file does not contain a column named 'x'. Columns found: ",
       paste(colnames(primate_df), collapse = ", "))
}
primate_cluster_2 <- unique(stats::na.omit(primate_df$x))

# -----------------------------
# Download + prepare GEO datasets
# -----------------------------
message("Downloading GEO datasets...")
gse30660 <- getGEO("GSE30660", GSEMatrix = TRUE)
gse82137 <- getGEO("GSE82137", GSEMatrix = TRUE)

expr_30660 <- exprs(gse30660[[1]])
expr_82137 <- exprs(gse82137[[1]])

# Match your original behaviour: drop columns 4:15 for GSE82137 if present
if (ncol(expr_82137) >= max(gse82137_drop_cols)) {
  expr_82137 <- expr_82137[, -gse82137_drop_cols, drop = FALSE]
}

# Phenotype data (used for conditions)
pheno_30660 <- pData(gse30660[[1]])
pheno_82137 <- pData(gse82137[[1]])

# -----------------------------
# Convert probes -> HGNC symbols
# -----------------------------
message("Mapping probes to HGNC symbols and collapsing duplicates...")
expr_30660_sym <- probe_to_symbol_matrix(expr_30660)
expr_82137_sym <- probe_to_symbol_matrix(expr_82137)

# -----------------------------
# Load NORM bulk expression and convert Ensembl -> HGNC
# -----------------------------
# Your original code used read.csv with sep='' which is fragile.
# We'll use read.delim with check.names=FALSE to preserve sample names.
message("Loading NORM expression matrix...")
norm_df <- read.delim(
  norm_expr_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# If the file contains a first column with gene IDs instead of rownames, handle it.
# Common patterns: first column named "gene" or "Gene" or blank.
if (!is.null(norm_df[[1]]) && !all(grepl("^[0-9.]+$", norm_df[[1]]))) {
  # If first column looks like IDs, set rownames from it.
  rownames(norm_df) <- norm_df[[1]]
  norm_df <- norm_df[, -1, drop = FALSE]
}

# Ensure numeric matrix
norm_df[] <- lapply(norm_df, function(x) as.numeric(as.character(x)))
if (anyNA(norm_df)) {
  warning("NORM expression contains NA after numeric conversion. Check input file formatting.")
}

message("Converting Ensembl IDs to HGNC symbols (biomaRt)...")
expr_norm_sym <- ensembl_to_hgnc_matrix(norm_df)

# -----------------------------
# Run GSVA
# -----------------------------
message("Running GSVA...")
gsva_30660 <- run_gsva_single_set(expr_30660_sym, primate_cluster_2, verbose = TRUE)
gsva_82137 <- run_gsva_single_set(expr_82137_sym, primate_cluster_2, verbose = TRUE)

# For RNA-seq style continuous expression, Gaussian kcdf is often used
gsva_norm <- run_gsva_single_set(expr_norm_sym, primate_cluster_2, kcdf = "Gaussian", verbose = TRUE)

# -----------------------------
# Plot: GSE30660
# -----------------------------
# Your original code used: pheno_ALI_1$`Smoker status:ch1`
# We'll try to find that column, otherwise fail with a helpful message.
col_30660 <- "Smoker status:ch1"
if (!col_30660 %in% colnames(pheno_30660)) {
  stop(
    "Expected phenotype column not found in GSE30660: `", col_30660, "`\n",
    "Available columns include: ", paste(colnames(pheno_30660), collapse = ", "), "\n",
    "Update `col_30660` to the correct column name for smoker status."
  )
}

cond_30660 <- pheno_30660[[col_30660]]
p_30660 <- make_gsva_boxplot(
  gsva_vec = as.numeric(gsva_30660[1, ]),
  condition = cond_30660,
  ylab = "GSVA score (primate_cluster_2)",
  title = "GSE30660"
)
print(p_30660)

if (isTRUE(save_plots)) {
  save_plot(p_30660, "cluster_2_in_GSE30660.tiff", width = 4, height = 3)
}

# -----------------------------
# Plot: GSE82137
# -----------------------------
# Your original code used: pheno_ALI_2$smoke
col_82137 <- "smoke"
if (!col_82137 %in% colnames(pheno_82137)) {
  stop(
    "Expected phenotype column not found in GSE82137: `", col_82137, "`\n",
    "Available columns include: ", paste(colnames(pheno_82137), collapse = ", "), "\n",
    "Update `col_82137` to the correct column name for smoke status."
  )
}

cond_82137 <- pheno_82137[[col_82137]]
p_82137 <- make_gsva_boxplot(
  gsva_vec = as.numeric(gsva_82137[1, ]),
  condition = cond_82137,
  ylab = "GSVA score (primate_cluster_2)",
  title = "GSE82137"
)
print(p_82137)

if (isTRUE(save_plots)) {
  save_plot(p_82137, "cluster_2_in_GSE82137.tiff", width = 4, height = 3)
}

# -----------------------------
# Plot: NORM (with stats)
# -----------------------------
# Your original code assumes `sample_conditions_NORM` exists.
# For GitHub reproducibility, we load it from a file if present,
# otherwise stop with instructions.
if (file.exists(norm_conditions_file)) {
  cond_df <- read.csv(norm_conditions_file, stringsAsFactors = FALSE)
  if (!all(c("sample", "Condition") %in% colnames(cond_df))) {
    stop("NORM conditions file must contain columns: sample, Condition")
  }
  
  # Align conditions to expression columns
  samples <- colnames(expr_norm_sym)
  cond_df <- cond_df %>% filter(sample %in% samples)
  
  if (nrow(cond_df) != length(samples)) {
    stop(
      "NORM conditions file does not cover all samples.\n",
      "Expected ", length(samples), " samples; found conditions for ", nrow(cond_df), "."
    )
  }
  
  # Reorder to match expression columns
  cond_df <- cond_df[match(samples, cond_df$sample), , drop = FALSE]
  cond_norm <- factor(cond_df$Condition, levels = c("Smoker", "Nonsmoker"))
} else {
  stop(
    "Missing NORM conditions file: ", norm_conditions_file, "\n",
    "Create it to make the script reproducible.\n",
    "Format: CSV with columns `sample` and `Condition` (Smoker/Nonsmoker)."
  )
}

gsva_scores_norm <- data.frame(
  GSVA_Score = as.numeric(gsva_norm[1, ]),
  Condition = cond_norm
)

p_norm <- ggplot(gsva_scores_norm, aes(x = Condition, y = GSVA_Score, fill = Condition)) +
  geom_boxplot(color = "black", width = 0.6) +
  geom_jitter(width = 0.2, size = 2, color = "black") +
  labs(
    title = NULL,
    x = NULL,
    y = "GSVA score of human lung evolution",
    fill = NULL
  ) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  stat_compare_means(
    comparisons = list(c("Smoker", "Nonsmoker")),
    method = "t.test",
    label = "p.format"
  )

print(p_norm)

if (isTRUE(save_plots)) {
  ggsave(
    filename = file.path(plot_dir, "cluster_2_in_NORM.tiff"),
    plot = p_norm,
    width = 4,
    height = 5,
    units = "in",
    dpi = plot_dpi,
    bg = plot_bg
  )
}

message("Done. Figures saved to: ", plot_dir)
