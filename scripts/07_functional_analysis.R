## 07_functional_analysis.R
## Functional annotation and enrichment analysis of DE genes
## Uses clusterProfiler for GO and KEGG ORA
## Follows the approach from BINF*6110 Lecture 8 tutorial

# ---- Load packages ----
library(DESeq2)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(ggplot2)
library(dplyr)

# ---- Load workspace from DESeq2 analysis ----
base_dir <- "C:/Users/yazan/Desktop/BINF6110-Assignment2"
output_dir <- file.path(base_dir, "output_files", "R_analysis")
load(file.path(output_dir, "deseq2_workspace.RData"))

# ---- Gene ID mapping ----
# Our gene IDs are already yeast ORF names (e.g. YAL001C)
# Map to common gene names for labeling
all_gene_ids <- rownames(res_early_mature)

gene_map <- bitr(all_gene_ids,
                 fromType = "ORF",
                 toType = c("GENENAME", "ENTREZID"),
                 OrgDb = org.Sc.sgd.db)

cat("Mapped", nrow(gene_map), "out of", length(all_gene_ids), "genes\n")

# ---- Prepare gene lists for ORA (Early vs Mature comparison) ----
res_df <- as.data.frame(res_early_mature)
res_df$ORF <- rownames(res_df)
res_df <- merge(res_df, gene_map, by = "ORF", all.x = TRUE)

# Significant upregulated genes (higher in Mature vs Early)
up_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  pull(ORF) %>%
  unique()

# Significant downregulated genes (lower in Mature vs Early)
down_genes <- res_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  pull(ORF) %>%
  unique()

# All significant DE genes
sig_genes <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(ORF) %>%
  unique()

# Background gene set
all_genes <- res_df$ORF %>% unique()

cat("Upregulated genes:", length(up_genes), "\n")
cat("Downregulated genes:", length(down_genes), "\n")
cat("Total DE genes:", length(sig_genes), "\n")
cat("Background genes:", length(all_genes), "\n")

# ---- GO Enrichment ORA ----

# Biological Process
ego_bp <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   keyType = "ORF",
                   ont = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)

# Molecular Function
ego_mf <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   keyType = "ORF",
                   ont = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)

# Cellular Component
ego_cc <- enrichGO(gene = sig_genes,
                   universe = all_genes,
                   OrgDb = org.Sc.sgd.db,
                   keyType = "ORF",
                   ont = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = FALSE)

# ---- Separate ORA for up/downregulated (for compareCluster) ----
compare_go <- compareCluster(
  geneCluster = list(Upregulated = up_genes, Downregulated = down_genes),
  fun = "enrichGO",
  OrgDb = org.Sc.sgd.db,
  keyType = "ORF",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2
)

# ---- KEGG Enrichment ----
# KEGG needs ENTREZID for yeast (organism code: sce)
sig_entrez <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

all_entrez <- res_df %>%
  filter(!is.na(ENTREZID)) %>%
  pull(ENTREZID) %>%
  unique()

kegg_enrich <- enrichKEGG(gene = sig_entrez,
                          universe = all_entrez,
                          organism = "sce",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2)

# ---- GSEA (ranked gene list approach) ----
# Rank all genes by log2FoldChange
gsea_list <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(desc(log2FoldChange)) %>%
  distinct(ORF, .keep_all = TRUE)

gene_ranks <- setNames(gsea_list$log2FoldChange, gsea_list$ORF)
gene_ranks <- sort(gene_ranks, decreasing = TRUE)

gsea_go <- gseGO(geneList = gene_ranks,
                 OrgDb = org.Sc.sgd.db,
                 keyType = "ORF",
                 ont = "BP",
                 minGSSize = 15,
                 maxGSSize = 500,
                 pvalueCutoff = 0.05)

# ---- Export results ----
if (nrow(as.data.frame(ego_bp)) > 0)
  write.csv(as.data.frame(ego_bp), file.path(output_dir, "GO_BP_enrichment.csv"))
if (nrow(as.data.frame(ego_mf)) > 0)
  write.csv(as.data.frame(ego_mf), file.path(output_dir, "GO_MF_enrichment.csv"))
if (nrow(as.data.frame(ego_cc)) > 0)
  write.csv(as.data.frame(ego_cc), file.path(output_dir, "GO_CC_enrichment.csv"))
if (nrow(as.data.frame(kegg_enrich)) > 0)
  write.csv(as.data.frame(kegg_enrich), file.path(output_dir, "KEGG_enrichment.csv"))
if (nrow(as.data.frame(gsea_go)) > 0)
  write.csv(as.data.frame(gsea_go), file.path(output_dir, "GSEA_GO_BP.csv"))

# Save workspace for visualization script
save.image(file = file.path(output_dir, "functional_workspace.RData"))

cat("Functional analysis complete. Results saved to:", output_dir, "\n")
