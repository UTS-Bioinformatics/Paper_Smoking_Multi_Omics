# GSVA_smk_sig.R
# Purpose:
#   - Build a smoking signature gene set from your DE list
#   - Download multiple GEO datasets and map probe IDs -> HGNC symbols
#   - Run GSVA for the smoking signature across datasets
#   - Plot GSVA scores with statistical comparisons
#   - (Optional) Plot AhR ChIP-seq fold enrichment by simplified annotation
#
# Repo assumptions (relative to repo root):
#   - A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke/de_genelist_smoking_NORM_annotated.csv
#   - (Optional) data/chipseq/AhR_chipseq.csv
#   - Output figures saved to: results/figures/
#

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(GSVA)
  library(limma)
  library(GEOquery)
  library(hgu133plus2.db)
  library(illuminaHumanv4.db)
  library(AnnotationDbi)
  library(ggpubr)
})

# -----------------------------
# User-configurable parameters
# -----------------------------
geneset_dir <- "A_Comprehensive_Multi_omics_and_Functional_Study_of_Evolutionary_Adaptive-Responses_to_Smoke"
smoking_file <- file.path(geneset_dir, "de_genelist_smoking_NORM_annotated.csv")

# Optional ChIP-seq file (set to NULL to skip)
ahr_chip_file <- file.path("data", "chipseq", "AhR_chipseq.csv")

# Output figures
plot_dir <- file.path("results", "figures")
save_plots <- TRUE
plot_dpi <- 600
plot_bg <- "white"

# Smoking signature thresholds
padj_cutoff <- 0.05
log2fc_cutoff <- 1

# Dataset-specific column drops (to match your original scripts)
gse82137_drop_cols <- 4:15
gse109576_drop_cols <- 10:12

# IMPORTANT:
# Your original scripts hard-coded sample labels as vectors. For reproducibility,
# keep that behaviour but centralise it here so it's easy to verify/modify.
sample_info <- list(
  GSE30660 = c("Control","Control","Control","Control","Smoker","Smoker","Smoker","Smoker"),
  GSE82137 = c("Smoker","Smoker","Smoker","Control","Control","Control"),
  GSE38332 = c("Mock","Mock","Mock","siRNA","siRNA","siRNA","NRF2_depleted","NRF2_depleted","NRF2_depleted"),
  GSE109576 = c("DMSO","DMSO","DMSO","AhR_activated","AhR_activated","AhR_activated","AhR_inhibited","AhR_inhibited","AhR_inhibited")
)

# Comparisons for multi-group datasets
comparisons <- list(
  GSE30660 = list(c("Smoker", "Control")),
  GSE82137 = list(c("Smoker", "Control")),
  GSE38332 = list(c("Mock", "siRNA"), c("Mock", "NRF2_depleted"), c("siRNA", "NRF2_depleted")),
  GSE109576 = list(c("DMSO", "AhR_inhibited"), c("DMSO", "AhR_activated"), c("AhR_activated", "AhR_inhibited"))
)

# Plot ordering (factor levels)
level_orders <- list(
  GSE30660 = c("Smoker", "Control"),
  GSE82137 = c("Smoker", "Control"),
  GSE38332 = c("NRF2_depleted", "siRNA", "Mock"),
  GSE109576 = c("AhR_inhibited", "AhR_activated", "DMSO")
)

# NOTE on colours:
# For GitHub/paper portability, we avoid hard-coding colours (especially bright cyan/orange).
# If you want them, set use_manual_fill <- TRUE and provide a palette below.
use_manual_fill <- FALSE
manual_fill <- list(
  GSE30660 = c("Smoker" = "orange", "Control" = "cyan"),
  GSE82137 = c("Smoker" = "orange", "Control" = "cyan"),
  GSE38332 = c("NRF2_depleted" = "orange", "siRNA" = "skyblue", "Mock" = "grey70"),
  GSE109576 = c("AhR_inhibited" = "orange", "AhR_activated" = "green", "DMSO" = "cyan")
)

# -----------------------------
# Checks + setup
# -----------------------------
if (!dir.exists(geneset_dir)) {
  stop(
    "Directory not found: ", geneset_dir, "\n",
    "Tip: set your working directory to the GitHub repository root before running this script."
  )
}
if (!file.exists(smoking_file)) stop("Missing file: ", smoking_file)

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Helper functions
# -----------------------------
collapse_by_gene_mean <- function(expr_mat, gene_symbols) {
  # expr_mat: matrix, rows correspond to probe IDs (already subsetted)
  # gene_symbols: vector same length as nrow(expr_mat), values are HGNC symbols
  df <- as.data.frame(expr_mat)
  df$Gene <- gene_symbols
  
  df <- df %>%
    filter(!is.na(Gene), Gene != "") %>%
    group_by(Gene) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")
  
  out <- as.matrix(df[, -1, drop = FALSE])
  rownames(out) <- df$Gene
  out
}

map_gpl570_probes_to_symbols <- function(expr_mat) {
  # GPL570 = hgu133plus2
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
  
  collapse_by_gene_mean(expr_mat, symbols)
}

map_illumina_v4_probes_to_symbols <- function(expr_mat) {
  symbols <- mapIds(
    illuminaHumanv4.db,
    keys = rownames(expr_mat),
    column = "SYMBOL",
    keytype = "PROBEID",
    multiVals = "first"
  )
  
  keep <- !is.na(symbols)
  expr_mat <- expr_mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  
  collapse_by_gene_mean(expr_mat, symbols)
}

run_gsva <- function(expr_mat, gene_set, kcdf = "Gaussian", verbose = TRUE) {
  gene_set <- unique(stats::na.omit(gene_set))
  gs <- list(smk_sig = gene_set)
  
  # Filter gene set to present genes
  gs[[1]] <- gs[[1]][gs[[1]] %in% rownames(expr_mat)]
  if (length(gs[[1]]) == 0) stop("No smoking signature genes are present in this expression matrix.")
  
  GSVA::gsva(expr_mat, gs, method = "gsva", kcdf = kcdf, verbose = verbose)
}

fit_limma <- function(gsva_mat, condition_vec) {
  condition <- factor(condition_vec)
  design <- model.matrix(~ condition)
  fit <- lmFit(gsva_mat, design)
  fit <- eBayes(fit)
  topTable(fit, coef = 2, adjust = "BH")
}

make_boxplot <- function(gsva_scores, condition_vec, dataset_name, ylab = "GSVA score") {
  df <- data.frame(
    GSVA_Score = as.numeric(gsva_scores),
    Condition = factor(condition_vec, levels = level_orders[[dataset_name]])
  )
  
  p <- ggplot(df, aes(x = Condition, y = GSVA_Score, fill = Condition)) +
    geom_boxplot(color = "black", width = 0.6) +
    geom_jitter(width = 0.2, size = 2, color = "black") +
    labs(y = ylab, x = NULL, fill = NULL) +
    theme_classic()
  
  if (isTRUE(use_manual_fill) && dataset_name %in% names(manual_fill)) {
    p <- p + scale_fill_manual(values = manual_fill[[dataset_name]])
  }
  
  if (dataset_name %in% names(comparisons)) {
    p <- p + stat_compare_means(
      comparisons = comparisons[[dataset_name]],
      method = "t.test",
      p.adjust.method = if (length(comparisons[[dataset_name]]) > 1) "BH" else "none",
      label = "p.format",
      tip.length = 0.01
    )
  }
  
  p
}

save_plot <- function(p, filename, width = 4, height = 3.5) {
  if (!isTRUE(save_plots)) return(invisible(NULL))
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
# Build smoking signature gene set
# -----------------------------
smk_df <- read.csv(smoking_file, header = TRUE, stringsAsFactors = FALSE)

required_cols <- c("padj", "log2FoldChange", "hgnc_symbol")
missing_cols <- setdiff(required_cols, colnames(smk_df))
if (length(missing_cols) > 0) {
  stop("Smoking signature file missing required columns: ", paste(missing_cols, collapse = ", "))
}

smk_sig <- smk_df %>%
  filter(padj < padj_cutoff, log2FoldChange > log2fc_cutoff) %>%
  pull(hgnc_symbol) %>%
  na.omit() %>%
  unique()

if (length(smk_sig) == 0) stop("Smoking signature gene set is empty after filtering. Check thresholds/file.")

message("Smoking signature size after filtering: ", length(smk_sig), " genes")

# -----------------------------
# GSE30660 (GPL570)
# -----------------------------
message("Processing GSE30660...")
gse30660 <- getGEO("GSE30660", GSEMatrix = TRUE)
expr_30660 <- exprs(gse30660[[1]])

expr_30660_sym <- map_gpl570_probes_to_symbols(expr_30660)

if (ncol(expr_30660_sym) != length(sample_info$GSE30660)) {
  stop("Sample label length mismatch for GSE30660: expr has ", ncol(expr_30660_sym),
       " samples, but sample_info has ", length(sample_info$GSE30660))
}

gsva_30660 <- run_gsva(expr_30660_sym, smk_sig, kcdf = "Gaussian", verbose = TRUE)
results_30660 <- fit_limma(gsva_30660, sample_info$GSE30660)

p_30660 <- make_boxplot(
  gsva_scores = gsva_30660["smk_sig", ],
  condition_vec = sample_info$GSE30660,
  dataset_name = "GSE30660",
  ylab = "GSVA score (smoking signature)"
)

print(p_30660)
save_plot(p_30660, "smk_sig_in_GSE30660.tiff", width = 4, height = 3.5)

# -----------------------------
# GSE82137 (GPL570, drop cols 4:15)
# -----------------------------
message("Processing GSE82137...")
gse82137 <- getGEO("GSE82137", GSEMatrix = TRUE)
expr_82137 <- exprs(gse82137[[1]])

if (ncol(expr_82137) >= max(gse82137_drop_cols)) {
  expr_82137 <- expr_82137[, -gse82137_drop_cols, drop = FALSE]
}

expr_82137_sym <- map_gpl570_probes_to_symbols(expr_82137)

if (ncol(expr_82137_sym) != length(sample_info$GSE82137)) {
  stop("Sample label length mismatch for GSE82137: expr has ", ncol(expr_82137_sym),
       " samples, but sample_info has ", length(sample_info$GSE82137))
}

gsva_82137 <- run_gsva(expr_82137_sym, smk_sig, kcdf = "Gaussian", verbose = TRUE)
results_82137 <- fit_limma(gsva_82137, sample_info$GSE82137)

p_82137 <- make_boxplot(
  gsva_scores = gsva_82137["smk_sig", ],
  condition_vec = sample_info$GSE82137,
  dataset_name = "GSE82137",
  ylab = "GSVA score (smoking signature)"
)

print(p_82137)
save_plot(p_82137, "smk_sig_in_GSE82137.tiff", width = 4, height = 3.5)

# -----------------------------
# GSE38332 (GPL570 via hgu133plus2; your original used GPL570 annotation table)
# -----------------------------
message("Processing GSE38332...")
gse38332 <- getGEO("GSE38332", GSEMatrix = TRUE)
expr_38332 <- exprs(gse38332[[1]])

# Safer + faster than pulling full GPL table each time:
expr_38332_sym <- map_gpl570_probes_to_symbols(expr_38332)

if (ncol(expr_38332_sym) != length(sample_info$GSE38332)) {
  stop("Sample label length mismatch for GSE38332: expr has ", ncol(expr_38332_sym),
       " samples, but sample_info has ", length(sample_info$GSE38332))
}

gsva_38332 <- run_gsva(expr_38332_sym, smk_sig, kcdf = "Gaussian", verbose = TRUE)
results_38332 <- fit_limma(gsva_38332, sample_info$GSE38332)

p_38332 <- make_boxplot(
  gsva_scores = gsva_38332["smk_sig", ],
  condition_vec = sample_info$GSE38332,
  dataset_name = "GSE38332",
  ylab = "GSVA score (smoking signature)"
)

print(p_38332)
save_plot(p_38332, "smk_sig_in_GSE38332.tiff", width = 5, height = 3.5)

# -----------------------------
# GSE109576 (Illumina HumanHT-12 v4; drop cols 10:12)
# -----------------------------
message("Processing GSE109576...")
gse109576 <- getGEO("GSE109576", GSEMatrix = TRUE)
expr_109576 <- exprs(gse109576[[1]])

if (ncol(expr_109576) >= max(gse109576_drop_cols)) {
  expr_109576 <- expr_109576[, -gse109576_drop_cols, drop = FALSE]
}

expr_109576_sym <- map_illumina_v4_probes_to_symbols(expr_109576)

if (ncol(expr_109576_sym) != length(sample_info$GSE109576)) {
  stop("Sample label length mismatch for GSE109576: expr has ", ncol(expr_109576_sym),
       " samples, but sample_info has ", length(sample_info$GSE109576))
}

gsva_109576 <- run_gsva(expr_109576_sym, smk_sig, kcdf = "Gaussian", verbose = TRUE)
results_109576 <- fit_limma(gsva_109576, sample_info$GSE109576)

p_109576 <- make_boxplot(
  gsva_scores = gsva_109576["smk_sig", ],
  condition_vec = sample_info$GSE109576,
  dataset_name = "GSE109576",
  ylab = "GSVA score (smoking signature)"
)

print(p_109576)
save_plot(p_109576, "smk_sig_in_GSE109576.tiff", width = 5, height = 3.5)

# -----------------------------
# Optional: AhR ChIP-seq plot
# -----------------------------
if (!is.null(ahr_chip_file) && file.exists(ahr_chip_file)) {
  message("Processing AhR ChIP-seq: ", ahr_chip_file)
  
  chip <- read.csv(ahr_chip_file, stringsAsFactors = FALSE)
  
  required_cols <- c("Annotation", "fold_enrichment")
  missing <- setdiff(required_cols, colnames(chip))
  if (length(missing) > 0) {
    stop("AhR chip-seq file missing required columns: ", paste(missing, collapse = ", "))
  }
  
  chip <- chip %>%
    mutate(
      Simplified_Annotation = case_when(
        grepl("Intergenic", Annotation, ignore.case = TRUE) ~ "Intergenic",
        grepl("promoter-TSS", Annotation, ignore.case = TRUE) ~ "Promoter_TSS",
        grepl("intron", Annotation, ignore.case = TRUE) ~ "Intron",
        TRUE ~ "Other"
      )
    )
  
  chip$Simplified_Annotation <- factor(
    chip$Simplified_Annotation,
    levels = c("Intergenic", "Intron", "Promoter_TSS", "Other")
  )
  
  chip_comparisons <- list(
    c("Intergenic", "Intron"),
    c("Promoter_TSS", "Intron"),
    c("Intergenic", "Promoter_TSS")
  )
  
  p_chip <- ggplot(chip, aes(x = Simplified_Annotation, y = fold_enrichment, fill = Simplified_Annotation)) +
    geom_boxplot(color = "black", width = 0.6) +
    geom_jitter(width = 0.2, size = 2, color = "black") +
    labs(x = NULL, y = "Fold enrichment") +
    theme_classic() +
    theme(legend.position = "none") +
    stat_compare_means(
      comparisons = chip_comparisons,
      method = "t.test",
      p.adjust.method = "BH",
      label = "p.format",
      tip.length = 0.01
    )
  
  print(p_chip)
  save_plot(p_chip, "AhR_chipseq.tiff", width = 5, height = 3.5)
  
} else {
  message("Skipping AhR ChIP-seq plot (file not found). If you want it, add: ", ahr_chip_file)
}

# -----------------------------
# Optional: write limma results tables
# -----------------------------
# These are helpful for paper supplements / reproducibility.
write.csv(results_30660, file.path("results", "tables_GSE30660_smk_sig_limma.csv"), row.names = TRUE)
write.csv(results_82137, file.path("results", "tables_GSE82137_smk_sig_limma.csv"), row.names = TRUE)
write.csv(results_38332, file.path("results", "tables_GSE38332_smk_sig_limma.csv"), row.names = TRUE)
write.csv(results_109576, file.path("results", "tables_GSE109576_smk_sig_limma.csv"), row.names = TRUE)

message("Done. Figures saved to: ", plot_dir)