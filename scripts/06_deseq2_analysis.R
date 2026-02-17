## 06_deseq2_analysis.R
## Differential expression analysis of yeast biofilm RNA-seq data
## Dataset: Mardanov et al. (2020) - S. cerevisiae flor yeast velum development
## Uses tximport for Salmon import and DESeq2 for DE analysis

# ---- Load packages ----
library(DESeq2)
library(tximport)
library(GenomicFeatures)
library(dplyr)

# ---- Set paths ----
base_dir <- "C:/Users/yazan/Desktop/BINF6110-Assignment2"
quant_dir <- file.path(base_dir, "salmon_quant")

# ---- Sample metadata ----
sample_table <- data.frame(
  sample_id = c("IL20", "IL21", "IL22", "IL23", "IL24", "IL25", "IL29", "IL30", "IL31"),
  run = c("SRR10551665", "SRR10551664", "SRR10551663",
          "SRR10551662", "SRR10551661", "SRR10551660",
          "SRR10551659", "SRR10551658", "SRR10551657"),
  stage = factor(c(rep("Early", 3), rep("Thin", 3), rep("Mature", 3)),
                 levels = c("Early", "Thin", "Mature")),
  days = c(rep(38, 3), rep(83, 3), rep(109, 3))
)
rownames(sample_table) <- sample_table$sample_id

# ---- Build tx2gene mapping from GTF ----
gtf_file <- file.path(base_dir, "GCF_000146045.2_R64_genomic.gtf")
txdb <- makeTxDbFromGFF(gtf_file)
k <- keys(txdb, keytype = "TXNAME")
tx2gene <- AnnotationDbi::select(txdb, k, "GENEID", "TXNAME")

# ---- Import Salmon counts with tximport ----
files <- file.path(quant_dir, sample_table$run, "quant.sf")
names(files) <- sample_table$sample_id

# Verify all files exist
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

# Quick look at the count matrix
cat("Dimensions of count matrix:", dim(txi$counts), "\n")
cat("First few genes:\n")
head(txi$counts)

# ---- Create DESeq2 dataset ----
dds <- DESeqDataSetFromTximport(txi, colData = sample_table, design = ~ stage)

# Pre-filter low count genes (at least 10 reads total across all samples)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
cat("Genes after filtering:", nrow(dds), "\n")

# ---- Run DESeq2 ----
dds <- DESeq(dds)
resultsNames(dds)

# ---- Pairwise comparisons ----
# Early vs Thin
res_early_thin <- results(dds, contrast = c("stage", "Thin", "Early"))
res_early_thin <- lfcShrink(dds, contrast = c("stage", "Thin", "Early"), type = "ashr")
summary(res_early_thin)

# Thin vs Mature
res_thin_mature <- results(dds, contrast = c("stage", "Mature", "Thin"))
res_thin_mature <- lfcShrink(dds, contrast = c("stage", "Mature", "Thin"), type = "ashr")
summary(res_thin_mature)

# Early vs Mature
res_early_mature <- results(dds, contrast = c("stage", "Mature", "Early"))
res_early_mature <- lfcShrink(dds, contrast = c("stage", "Mature", "Early"), type = "ashr")
summary(res_early_mature)

# ---- Likelihood Ratio Test (time-course approach) ----
# Tests for any gene expression change across all stages
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)
res_lrt <- results(dds_lrt)
summary(res_lrt)

cat("LRT significant genes (padj < 0.05):", sum(res_lrt$padj < 0.05, na.rm = TRUE), "\n")

# ---- Export DE results ----
output_dir <- file.path(base_dir, "output_files", "R_analysis")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(as.data.frame(res_early_thin), file.path(output_dir, "DE_early_vs_thin.csv"))
write.csv(as.data.frame(res_thin_mature), file.path(output_dir, "DE_thin_vs_mature.csv"))
write.csv(as.data.frame(res_early_mature), file.path(output_dir, "DE_early_vs_mature.csv"))
write.csv(as.data.frame(res_lrt), file.path(output_dir, "DE_LRT_all_stages.csv"))

# ---- Save workspace for downstream scripts ----
save(dds, dds_lrt, txi, sample_table, tx2gene,
     res_early_thin, res_thin_mature, res_early_mature, res_lrt,
     file = file.path(output_dir, "deseq2_workspace.RData"))

cat("DESeq2 analysis complete. Results saved to:", output_dir, "\n")
