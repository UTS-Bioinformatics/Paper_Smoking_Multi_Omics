# smk_sig_DGE.R
# Purpose:
#   - Load NORM expression matrix (exp) and metadata (info)
#   - Run edgeR QLF differential expression:
#       current_smoker vs never_smoker
#     adjusting for age and gender
#
# Data provenance:
#   - exp/info were extracted from the public GEO dataset: GSE237252
#   - Upstream steps (e.g., gene ID conversion, filtering) have already been applied
#     and are reflected in the provided exp/info files.
#
# Repo assumptions (relative to repo root):
#   - NORM_exp.csv   (rows=genes, cols=samples)
#   - NORM_meta.csv  (rows=samples; must include smoking_status, age, gender)
#
# Output:
#   - results/tables/baseline.csv

suppressPackageStartupMessages({
  library(edgeR)
})

# -----------------------------
# Load exp and info
# -----------------------------
exp_file  <- file.path("NORM_exp.csv")
info_file <- file.path("NORM_meta.csv")

if (!file.exists(exp_file))  stop("Missing file: ", exp_file)
if (!file.exists(info_file)) stop("Missing file: ", info_file)

# Expression matrix: rows=genes, cols=samples
exp <- read.csv(exp_file, row.names = 1, check.names = FALSE)

# Metadata: rows=samples (rnaseq.id or sample IDs)
info <- read.csv(info_file, row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)

# Minimal required columns
required_cols <- c("smoking_status", "age", "gender")
missing_cols <- setdiff(required_cols, colnames(info))
if (length(missing_cols) > 0) {
  stop("Metadata file missing required columns: ", paste(missing_cols, collapse = ", "))
}

# Ensure types
info$smoking_status <- factor(info$smoking_status)
info$age <- as.numeric(info$age)
info$gender <- factor(info$gender)

# Align samples between exp and info (keep common, in same order)
common_samples <- intersect(colnames(exp), rownames(info))
if (length(common_samples) == 0) stop("No overlapping sample IDs between exp columns and info rownames.")

exp <- exp[, common_samples, drop = FALSE]
info <- info[common_samples, , drop = FALSE]

# Ensure exp is numeric matrix
exp <- as.matrix(exp)
storage.mode(exp) <- "numeric"

# -----------------------------
# Run edgeR process #
# -----------------------------
DGE <- DGEList(exp, samples = rownames(info), group = info$smoking_status)
keep.exprs <- filterByExpr(DGE)
DGE_after <- DGE[keep.exprs, keep.lib.sizes = F]
DGE_after <- calcNormFactors(DGE_after, method = "TMM")
design <- model.matrix(~0 + smoking_status + age + gender, data = info)
colnames(design)
DGE_after <- estimateDisp(DGE_after, design)
fit <- glmQLFit(DGE_after, design)
cont.matrix <- makeContrasts(smoking_statuscurrent_smoker - smoking_statusnever_smoker,
                             levels= design)
qlf <- glmQLFTest(fit, contrast = cont.matrix)
tT <- topTags(qlf, n = nrow(DGE_after))$table
tT$gene <- rownames(tT)

# Save to repo-friendly path
out_dir <- file.path("results", "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(tT, file.path(out_dir, "baseline.csv"), row.names = FALSE)
