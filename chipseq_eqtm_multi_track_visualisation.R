# chipseq_eqtm_multi_track_visualisation.R
# Purpose:
#   - Integrate ChIP-seq (AHR, AHRR, NRF2), DNA methylation (eQTM),
#     and gene annotation tracks
#   - Generate locus-level multi-track visualisations for each gene
#
# Data provenance:
#   - All ChIP-seq, methylation, and eQTM tables were derived from
#     public datasets described in the smoking-related study.
#   - These processed tables are included in this repository for
#     visualisation purposes only.
#.  - All data are provided in chiseq zip file.
#
# IMPORTANT NOTE ON DEPENDENCIES:
#   - The `GenomeGraphs` package is ONLY compatible with R < 4.0.0
#   - If you are using R >= 4.0.0, please create a separate R < 4.0.0
#     environment (e.g. via renv or conda) to run this script.
#
# Repo assumptions (relative to repo root):
#   data/chipseq/
#     ├─ AHR_24hr.txt
#     ├─ AHR_45min.txt
#     ├─ AHRR_24hr.txt
#     ├─ NRF2_24hr.txt
#     ├─ gene_positions.txt
#     ├─ methylation_locations.txt
#     ├─ gene_id_to_symbol.txt
#     ├─ eqtm_gene_cpg_table.txt
#     └─ eqtm_all_significance_table.txt
#
# Output:
#   - result_plots/<GENE_SYMBOL>.pdf   (one file per gene)
#
# How to run:
#   - Use R < 4.0.0

# -----------------------------
# Package loading
# -----------------------------
suppressPackageStartupMessages({
  library(ggfortify)
  library(GenomeGraphs)  # requires R < 4.0.0
  library(cowplot)
  library(gridExtra)
  library(ggplot2)
  library(ggrepel)
  library(biomaRt)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  library(BSgenome.Hsapiens.UCSC.hg19)
  library(Gviz)
  library(ggbio)
  library(rtracklayer)
  library(GenomicRanges)
  library(Homo.sapiens)
})

# -----------------------------
# Paths
# -----------------------------
data_dir <- file.path("data", "chipseq")
plot_dir <- "result_plots"

dir.create(plot_dir, showWarnings = FALSE)

# -----------------------------
# Load input files
# -----------------------------
AHR24file        <- read.delim(file.path(data_dir, "AHR_24hr.txt"))
AHR45minfile     <- read.delim(file.path(data_dir, "AHR_45min.txt"))
AHRR24hrs_file   <- read.delim(file.path(data_dir, "AHRR_24hr.txt"))
NRF2_file        <- read.delim(file.path(data_dir, "NRF2_24hr.txt"))
genepos          <- read.delim(file.path(data_dir, "gene_positions.txt"))
methyl_locations <- read.delim(file.path(data_dir, "methylation_locations.txt"))
gene_symbols     <- read.delim(file.path(data_dir, "gene_id_to_symbol.txt"))
eQTMfile         <- read.delim(file.path(data_dir, "eqtm_gene_cpg_table.txt"), sep = " ")
all_eQTMfile     <- read.delim(file.path(data_dir, "eqtm_all_significance_table.txt"), sep = " ")

chipseq_file <- read.delim(file.path(data_dir, "final_chipseq_table.txt"))
chipseq_file <- chipseq_file[
  !(is.na(chipseq_file$start) &
      is.na(chipseq_file$start.1) &
      is.na(chipseq_file$start.2)),
]

rownames(methyl_locations) <- methyl_locations[, 1]

# -----------------------------
# Prepare gene list
# -----------------------------
chipseq_genes <- unique(as.character(chipseq_file$gene))
chipseq_genes <- chipseq_genes[chipseq_genes != ""]

# -----------------------------
# Loop over genes
# -----------------------------
for (i in seq_along(chipseq_genes)) {
  
  gene_id <- chipseq_genes[i]
  genename <- unique(gene_symbols$GENENAME[gene_symbols$ID == gene_id])
  
  if (length(genename) == 0 || genename == "") next
  
  # Gene position
  gene_pos <- genepos[genepos$ID == gene_id, ]
  chr <- as.character(gene_pos$chr[1])
  
  # Extract eQTM CpGs
  cpg_gene <- all_eQTMfile[as.character(all_eQTMfile[, 2]) == gene_id, ]
  cpg_gene <- cbind(cpg_gene, methyl_locations[as.character(cpg_gene[, 1]), ])
  
  # Significance annotation
  cpg_gene$Significant <- ifelse(
    as.numeric(cpg_gene$FDR) < 0.05 & as.numeric(cpg_gene$statistic) < 0, "Decreased",
    ifelse(as.numeric(cpg_gene$FDR) < 0.05 & as.numeric(cpg_gene$statistic) > 0, "Increased",
           "Not Sig")
  )
  
  # Genomic window
  center <- as.numeric(gene_pos$pos[1])
  flank <- 5000
  start <- center - flank
  end   <- center + flank
  
  # -----------------------------
  # Plot 1: eQTM Manhattan-style
  # -----------------------------
  p1 <- ggplot(cpg_gene, aes(x = as.numeric(pos), y = -log10(as.numeric(pvalue)))) +
    geom_point(aes(color = Significant)) +
    geom_hline(yintercept = -log10(0.0008), linetype = "dashed") +
    geom_text_repel(
      data = subset(cpg_gene, as.numeric(FDR) < min(as.numeric(FDR)) + 0.0005),
      aes(label = snps),
      size = 3
    ) +
    xlim(c(start, end)) +
    theme_bw() +
    theme(legend.position = "none")
  
  # -----------------------------
  # Plot 2: ChIP-seq tracks
  # -----------------------------
  extract_chip <- function(df, label) {
    df[df$chr == chr & df$start > start & df$start < end, ] |>
      transform(type = label)
  }
  
  chip_all <- rbind(
    extract_chip(AHR24file, "AHR_24hr"),
    extract_chip(AHR45minfile, "AHR_45min"),
    extract_chip(NRF2_file, "NRF2_24hr"),
    extract_chip(AHRR24hrs_file, "AHRR_24hr")
  )
  
  p2 <- ggplot(chip_all) +
    geom_segment(aes(x = start, xend = end, y = type, yend = type,
                     color = fold_enrichment), size = 3) +
    xlim(c(start, end)) +
    theme_void()
  
  # -----------------------------
  # Plot 3: Gene annotation
  # -----------------------------
  data(genesymbol, package = "biovizBase")
  gr <- GRanges(chr, IRanges(start, end))
  gene_track <- range(genesymbol[
    genesymbol$ensembl_id == gene_id |
      genesymbol$symbol == genename
  ])
  
  pdf(file.path(plot_dir, paste0(genename, ".pdf")), width = 10, height = 6)
  
  p_tx <- autoplot(Homo.sapiens, which = gene_track)
  print(tracks(
    Gene = p_tx,
    ChipSeq = p2,
    eQTM = p1,
    heights = c(3)
  ))
  
  dev.off()
  
  message("Processed gene: ", genename, " (", i, "/", length(chipseq_genes), ")")
}

# This script / directory contains a vendored copy of the R package **GenomeGraphs**,
# which was originally distributed via Bioconductor.

# GenomeGraphs is no longer available through standard Bioconductor repositories.
# It is included here solely to ensure reproducibility of analyses in the
# associated manuscript.

#   - Original authorship and license are preserved
#   - No claim of authorship or ownership is made
#   - The package is used without modification

# If possible, users should prefer an official Bioconductor source.
