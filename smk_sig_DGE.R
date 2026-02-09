library(biomaRt)
library(dplyr)
library(tidyr)
library(edgeR)

# load NORM data #
gene <- read.table("/Users/XXX/1504_VanenBerge.expression.genelevel.v75.htseq.txt", 
                   header = T, row.names = 1)
info <- read.csv("/Users/XXX/Database_biopten_v2_6.csv")
info$gender <- factor(info$gender)
info$age <- as.numeric(info$age)

# convert to gene symbols #
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- rownames(gene)
G_list <- getBM(filters= "ensembl_gene_id", 
                attributes= c("ensembl_gene_id","hgnc_symbol"), 
                values = genes, mart = mart)
genes <- as.matrix(rownames(gene))
data <- merge(genes, G_list, by.x = "V1", by.y = "ensembl_gene_id")
exp <- as.matrix(gene[data[,1],])
rownames(exp) <- data[,2]

# remove asthma patients #
info <- info[which(info$asthma == "healthy"),]
rownames(info) <- info$rnaseq.id
exp <- exp[, colnames(exp) %in% rownames(info), drop = T]

any(duplicated(rownames(exp)))
which(duplicated(rownames(exp)))
exp <- exp[!duplicated(rownames(exp)), , drop = FALSE]

# handling clinical data #
info <- info %>%
  separate(
    col   = normgroup,
    into  = c("age_group", "smoking_status"),
    sep   = "\\.",        
    remove = FALSE        
  )
info <- info[rownames(info) %in% colnames(exp), , drop = T]
exp <- exp[, rownames(info)]

# QC #
libsize <- colSums(exp)
n_expressed <- colSums(exp > 1)
gene_means <- rowMeans(exp)
gene_vars  <- apply(exp, 1, var)
log_exp <- log2(exp + 1)
sample_dist <- dist(t(log_exp))
hc <- hclust(sample_dist)

QC <- c(
  "X102_NORM_",
  "X104_NORM_",
  "X105_NORM_",
  "X108_NORM_",
  "X109_NORM_",
  "X10_NORM_",
  "X111_NORM_",
  "X113_NORM_",
  "X114_NORM_",
  "X122_NORM_",
  "X123_NORM_",
  "X124_NORM_",
  "X125_NORM_",
  "X128_NORM_",
  "X129_NORM_",
  "X130_NORM_",
  "X135_NORM_",
  "X136_NORM_",
  "X137_NORM_",
  "X138_NORM_",
  "X139_NORM_",
  "X145_NORM_",
  "X146_NORM_",
  "X147_NORM_",
  "X148_NORM_",
  "X149_NORM_",
  "X153_NORM_",
  "X154_NORM_",
  "X155_NORM_",
  "X158_NORM_",
  "X159_NORM_",
  "X15_NORM_",
  "X161_NORM_",
  "X163_NORM_",
  "X164_NORM_",
  "X166_NORM_",
  "X16_NORM_",
  "X17_NORM_",
  "X18_NORM_",
  "X20_NORM_",
  "X22_NORM_",
  "X24_NORM_",
  "X25_NORM_",
  "X28_NORM_",
  "X2_NORM_",
  "X31_NORM_",
  "X32_NORM_",
  "X33_NORM_",
  "X37_NORM_",
  "X39_NORM_",
  "X40_NORM_",
  "X41_NORM_",
  "X42_NORM_",
  "X46_NORM_",
  "X47_NORM_",
  "X48_NORM_",
  "X4_NORM_",
  "X55_NORM_",
  "X60_NORM_",
  "X63_NORM_",
  "X66_NORM_",
  "X68_NORM_",
  "X71_NORM_",
  "X74_NORM_",
  "X77_NORM_",
  "X78_NORM_",
  "X80_NORM_",
  "X81_NORM_",
  "X84_NORM_",
  "X85_NORM_",
  "X86_NORM_",
  "X88_NORM_",
  "X89_NORM_",
  "X8_NORM_",
  "X96_NORM_",
  "X98_NORM_",
  "X99_NORM_") # samples to keep

# need to remove (n=6): 23, 49, 106, 118, 52, 92 #
samples_to_remove <- c("X23_NORM_", "X49_NORM_", "X106_NORM_", 
                       "X118_NORM_", "X52_NORM_", "X92_NORM_")
exp <- exp[, !colnames(exp) %in% samples_to_remove]
info <- info[!rownames(info) %in% samples_to_remove,]

exp <- exp[, colnames(exp) %in% QC]
info <- info[rownames(info) %in% QC,]

sum(info$smoking_status == "non-smoker")
sum(info$smoking_status == "smoker")

info$smoking_status[info$smoking_status == "smoker"] <- "current-smoker"
info$smoking_status[info$smoking_status == "non-smoker"] <- "never-smoker"

info$smoking_status <- gsub("-", "_", info$smoking_status)
info$smoking_status <- factor(info$smoking_status)

# run edgeR process #
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
write.csv(tT, "/Users/XXX/smk_sig.csv", 
          row.names = FALSE)
