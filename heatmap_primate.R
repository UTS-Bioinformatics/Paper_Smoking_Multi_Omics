# heatmap_primate.R
# Purpose:
#   - Load primate (and never-smoker human) z-score matrix
#   - Standardise/rename columns (species + SRR runs -> Never_smoker_*)
#   - Draw a heatmap (pheatmap) and save as TIFF
#
# Repo assumptions (relative to repo root):
#   - data/Primate/primates_nonsmokers_zscores.txt
#   - Output saved to: results/figures/new_primate_heatmap.tiff
#

suppressPackageStartupMessages({
  library(pheatmap)
  library(viridis)
})

# -----------------------------
# User-configurable parameters
# -----------------------------
input_file <- file.path("data", "Primate", "primates_nonsmokers_zscores.txt")

plot_dir <- file.path("results", "figures")
output_file <- file.path(plot_dir, "new_primate_heatmap.tiff")

dpi <- 600
plot_width <- 8
plot_height <- 5
bg <- "white"

# Column naming assumptions:
#   - The first 11 columns correspond to the listed species (in that order)
#   - Additional columns with names starting with "SRR" are human never-smoker runs
species_names_11 <- c(
  "Baboon",
  "Chimpanzee_1",
  "Chimpanzee_2",
  "Cynomolgus_macaque_1",
  "Cynomolgus_macaque_2",
  "Marmoset",
  "Mouse_lemur",
  "Pigtailed_macaque",
  "Rhesus_macaque",
  "Sooty_mangabey",
  "Squirrel_monkey"
)

# -----------------------------
# Checks + setup
# -----------------------------
if (!file.exists(input_file)) {
  stop(
    "Missing input file: ", input_file, "\n",
    "Put it in the repo at that path or update `input_file`."
  )
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load data
# -----------------------------
dat <- read.delim(
  input_file,
  header = TRUE,
  row.names = 1,
  na.strings = c("N/A"),
  check.names = FALSE
)

# -----------------------------
# Rename columns
# -----------------------------
# 1) Rename SRR* columns to Never_smoker_#
srr_idx <- grep("^SRR", colnames(dat))
if (length(srr_idx) > 0) {
  colnames(dat)[srr_idx] <- paste0("Never_smoker_", seq_along(srr_idx))
} else {
  warning("No columns starting with 'SRR' found. Skipping Never_smoker_* renaming.")
}

# 2) Rename first 11 columns to species names (only if they exist)
if (ncol(dat) < length(species_names_11)) {
  warning(
    "Data has fewer than 11 columns (ncol=", ncol(dat), "). ",
    "Skipping species renaming for first 11 columns."
  )
} else {
  colnames(dat)[seq_along(species_names_11)] <- species_names_11
}

# -----------------------------
# Plot heatmap
# -----------------------------
mat <- as.matrix(dat)

# pheatmap returns a "pheatmap" object; ggsave() won't reliably save it.
# We save via tiff() + print().
breaks <- seq(-4, 4, length.out = 100)

p <- pheatmap(
  mat,
  cluster_cols = FALSE,
  fontsize_col = 10,
  fontsize_row = 4,
  na_col = "white",
  breaks = breaks,
  color = viridis(100),
  gaps_col = length(species_names_11),  # visual separation after the 11 species columns
  cutree_rows = 4,
  show_rownames = FALSE
)

# -----------------------------
# Save figure
# -----------------------------
tiff(
  filename = output_file,
  width = plot_width,
  height = plot_height,
  units = "in",
  res = dpi,
  compression = "lzw",
  bg = bg
)
print(p)
dev.off()

message("Saved heatmap to: ", output_file)