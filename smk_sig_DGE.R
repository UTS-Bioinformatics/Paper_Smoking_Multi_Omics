# smk_sig_DGE.R
# Purpose:
#   - Load NORM bronchial biopsy bulk RNA-seq gene-level counts and clinical metadata
#   - Convert Ensembl gene IDs to HGNC symbols
#   - Filter to healthy (non-asthma) participants and perform basic QC filtering
#   - Run edgeR quasi-likelihood differential expression:
#       current_smoker vs never_smoker
#     adjusting for age and gender
#   - Export a DE table that can be used to define a smoking gene signature
#
# Data provenance / availability:
#   - The NORM bulk RNA-seq dataset is publicly available under GEO accession:
#       GSE237252
#   - Therefore, the raw/count and clinical files are NOT provided in this GitHub repository.
#     Users should download the dataset from GEO (GSE237252) and place files locally, then
#     update the input paths below.
#
# Expected inputs (local paths; NOT in this repo):
#   - gene_count_file:  gene-level HTSeq count matrix (rows=Ensembl IDs, cols=samples)
#   - metadata_file:    clinical metadata table (must include rnaseq.id, asthma, normgroup, age, gender)
#
# Outputs:
#   - results/tables/de_genelist_smoking_NORM_annotated.csv
#

suppressPackageStartupMessages({
  library(biomaRt)
  library(dplyr)
  library(tidyr)
  library(edgeR)
})

# -----------------------------
# User-configurable parameters
# -----------------------------
# NOTE: These files are NOT included in the GitHub repo (data in GSE237252).
# Set these to your local file locations after downloading from GEO.
gene_count_file <- "PATH/TO/GSE237252_gene_counts.txt"   # e.g., gene-level HTSeq counts
metadata_file   <- "PATH/TO/GSE237252_clinical_metadata.csv"

# Output
out_dir  <- file.path("results", "tables")
out_file <- file.path(out_dir, "de_genelist_smoking_NORM_annotated.csv")

# Clinical filters
keep_asthma_value <- "healthy"   # keep only asthma == "healthy"

# QC: explicit allow-list of samples to keep (from your original script)
samples_keep_allowlist <- c(
  "X102_NORM_", "X104_NORM_", "X105_NORM_", "X108_NORM_", "X109_NORM_", "X10_NORM_",
  "X111_NORM_", "X113_NORM_", "X114_NORM_", "X122_NORM_", "X123_NORM_", "X124_NORM_",
  "X125_NORM_", "X128_NORM_", "X129_NORM_", "X130_NORM_", "X135_NORM_", "X136_NORM_",
  "X137_NORM_", "X138_NORM_", "X139_NORM_", "X145_NORM_", "X146_NORM_", "X147_NORM_",
  "X148_NORM_", "X149_NORM_", "X153_NORM_", "X154_NORM_", "X155_NORM_", "X158_NORM_",
  "X159_NORM_", "X15_NORM_",  "X161_NORM_", "X163_NORM_", "X164_NORM_", "X166_NORM_",
  "X16_NORM_",  "X17_NORM_",  "X18_NORM_",  "X20_NORM_",  "X22_NORM_",  "X24_NORM_",
  "X25_NORM_",  "X28_NORM_",  "X2_NORM_",   "X31_NORM_",  "X32_NORM_",  "X33_NORM_",
  "X37_NORM_",  "X39_NORM_",  "X40_NORM_",  "X41_NORM_",  "X42_NORM_",  "X46_NORM_",
  "X47_NORM_",  "X48_NORM_",  "X4_NORM_",   "X55_NORM_",  "X60_NORM_",  "X63_NORM_",
  "X66_NORM_",  "X68_NORM_",  "X71_NORM_",  "X74_NORM_",  "X77_NORM_",  "X78_NORM_",
  "X80_NORM_",  "X81_NORM_",  "X84_NORM_",  "X85_NORM_",  "X86_NORM_",  "X88_NORM_",
  "X89_NORM_",  "X8_NORM_",   "X96_NORM_",  "X98_NORM_",  "X99_NORM_"
)

# QC: explicit remove-list (from your original script)
samples_remove_list <- c("X23_NORM_", "X49_NORM_", "X106_NORM_", "X118_NORM_", "X52_NORM_", "X92_NORM_")

# -----------------------------
# Checks + setup
# -----------------------------
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(gene_count_file)) {
  stop(
    "Gene count file not found: ", gene_count_file, "\n",
    "This repo does not include the NORM bulk RNA-seq files. Download them from GSE237252 and update `gene_count_file`."
  )
}
if (!file.exists(metadata_file)) {
  stop(
    "Metadata file not found: ", metadata_file, "\n",
    "This repo does not include the NORM bulk RNA-seq files. Download them from GSE237252 and update `metadata_file`."
  )
}

# -----------------------------
# Load data
# -----------------------------
# Gene counts: rows = Ensembl gene IDs, columns = samples
gene_counts <- read.table(
  gene_count_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# Clinical metadata
info <- read.csv(metadata_file, stringsAsFactors = FALSE)

# Basic type fixes
if ("gender" %in% colnames(info)) info$gender <- factor(info$gender)
if ("age" %in% colnames(info)) info$age <- as.numeric(info$age)

required_meta_cols <- c("rnaseq.id", "asthma", "normgroup", "age", "gender")
missing_meta_cols <- setdiff(required_meta_cols, colnames(info))
if (length(missing_meta_cols) > 0) {
  stop("Metadata file missing required columns: ", paste(missing_meta_cols, collapse = ", "))
}

# -----------------------------
# Convert Ensembl -> HGNC symbols (biomaRt)
# -----------------------------
message("Converting Ensembl IDs to HGNC symbols (biomaRt)...")
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

ensembl_ids <- rownames(gene_counts)

gene_map <- getBM(
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  values = ensembl_ids,
  mart = mart
) %>%
  filter(!is.na(hgnc_symbol), hgnc_symbol != "")

if (nrow(gene_map) == 0) {
  stop("biomaRt returned no HGNC symbols. Check that rownames are Ensembl gene IDs and that you have internet access.")
}

# Join counts with mapping, then collapse duplicates by HGNC symbol (mean)
gene_counts$ensembl_gene_id <- rownames(gene_counts)
merged <- merge(gene_counts, gene_map, by = "ensembl_gene_id")

merged <- merged %>%
  select(-ensembl_gene_id) %>%
  group_by(hgnc_symbol) %>%
  summarise(across(where(is.numeric), mean), .groups = "drop")

exp <- as.matrix(merged[, -1, drop = FALSE])
rownames(exp) <- merged$hgnc_symbol

# -----------------------------
# Filter metadata: keep healthy only (remove asthma patients)
# -----------------------------
info <- info[info$asthma == keep_asthma_value, , drop = FALSE]
rownames(info) <- info$rnaseq.id

# Keep only samples present in both exp and info
common_samples <- intersect(colnames(exp), rownames(info))
exp <- exp[, common_samples, drop = FALSE]
info <- info[common_samples, , drop = FALSE]

# -----------------------------
# Parse normgroup into age_group + smoking_status
# -----------------------------
# normgroup looks like "<age_group>.<smoking_status>"
info <- info %>%
  separate(
    col = normgroup,
    into = c("age_group", "smoking_status"),
    sep = "\\.",
    remove = FALSE
  )

# -----------------------------
# QC filtering using your curated lists
# -----------------------------
# Remove explicit bad samples (if present)
exp <- exp[, !colnames(exp) %in% samples_remove_list, drop = FALSE]
info <- info[!rownames(info) %in% samples_remove_list, , drop = FALSE]

# Keep only curated QC allowlist (if present)
keep_samples <- intersect(colnames(exp), samples_keep_allowlist)
exp <- exp[, keep_samples, drop = FALSE]
info <- info[keep_samples, , drop = FALSE]

# -----------------------------
# Clean smoking status labels
# -----------------------------
# Original labels: "smoker" / "non-smoker"
info$smoking_status[info$smoking_status == "smoker"] <- "current_smoker"
info$smoking_status[info$smoking_status == "non-smoker"] <- "never_smoker"
info$smoking_status <- gsub("-", "_", info$smoking_status)
info$smoking_status <- factor(info$smoking_status)

message("N current_smoker: ", sum(info$smoking_status == "current_smoker", na.rm = TRUE))
message("N never_smoker:   ", sum(info$smoking_status == "never_smoker", na.rm = TRUE))

if (anyNA(info$smoking_status)) stop("smoking_status contains NA after parsing/recoding. Check `normgroup` formatting.")

# -----------------------------
# edgeR workflow (QLF)
# -----------------------------
message("Running edgeR...")

DGE <- DGEList(counts = exp, samples = info, group = info$smoking_status)

keep.exprs <- filterByExpr(DGE, group = info$smoking_status)
DGE <- DGE[keep.exprs, , keep.lib.sizes = FALSE]

DGE <- calcNormFactors(DGE, method = "TMM")

# Model: smoking_status + age + gender
design <- model.matrix(~ 0 + smoking_status + age + gender, data = info)

DGE <- estimateDisp(DGE, design)
fit <- glmQLFit(DGE, design)

# Contrast: current_smoker - never_smoker
cont <- makeContrasts(
  smoking_statuscurrent_smoker - smoking_statusnever_smoker,
  levels = design
)

qlf <- glmQLFTest(fit, contrast = cont)

tt <- topTags(qlf, n = nrow(DGE))$table
tt$gene <- rownames(tt)

# -----------------------------
# Save results
# -----------------------------
write.csv(tt, out_file, row.names = FALSE)
message("Saved DE table to: ", out_file)
