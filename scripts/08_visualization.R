## 08_visualization.R

# ---- Load packages ----
library(DESeq2)
library(clusterProfiler)
library(org.Sc.sgd.db)
library(enrichplot)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(dplyr)

# ---- Load workspaces ----
base_dir <- "C:/Users/yazan/Desktop/BINF6110-Assignment2"
output_dir <- file.path(base_dir, "output_files", "R_analysis")
load(file.path(output_dir, "functional_workspace.RData"))

fig_dir <- file.path(output_dir, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Figure 1: PCA Plot ----
vsd <- vst(dds, blind = FALSE)

pca_data <- plotPCA(vsd, intgroup = c("stage"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = stage, shape = stage)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("Early" = "#2C5F2D", "Thin" = "#E8A838", "Mature" = "#8B2252")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of Yeast Biofilm Samples by Development Stage") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figure1_pca.png"), pca_plot, width = 7, height = 6, dpi = 300)
cat("Figure 1 (PCA) saved.\n")

# ---- Figure 2: Volcano Plot (Early vs Mature) ----
res_df <- as.data.frame(res_early_mature)
res_df$ORF <- rownames(res_df)

# Map ORF to common gene names for labeling
gene_names <- bitr(res_df$ORF,
                   fromType = "ORF",
                   toType = "GENENAME",
                   OrgDb = org.Sc.sgd.db)
res_df <- merge(res_df, gene_names, by = "ORF", all.x = TRUE)

# Label genes as up/down/not significant
res_df$direction <- "Not Sig"
res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up"
res_df$direction[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"
res_df$direction <- factor(res_df$direction, levels = c("Up", "Down", "Not Sig"))

# Pick top genes to label
res_df_sig <- res_df %>%
  filter(direction != "Not Sig") %>%
  arrange(padj) %>%
  head(15)

volcano <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = direction)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "Not Sig" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  ggrepel::geom_text_repel(data = res_df_sig,
                           aes(label = ifelse(!is.na(GENENAME) & GENENAME != "", GENENAME, ORF)),
                           size = 3, max.overlaps = 20, color = "black") +
  labs(x = "Log2 Fold Change (Mature vs Early)",
       y = "-Log10 Adjusted P-value",
       title = "Differentially Expressed Genes: Mature vs Early Biofilm",
       color = "Direction") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "figure2_volcano.png"), volcano, width = 8, height = 7, dpi = 300)
cat("Figure 2 (Volcano) saved.\n")

# ---- Figure 3: Heatmap of Top DE Genes ----
# Top 30 DE genes from LRT (most variable across all stages)
res_lrt_ordered <- res_lrt[order(res_lrt$padj), ]
top_genes <- head(rownames(res_lrt_ordered), 30)

mat <- assay(vsd)[top_genes, ]

# Get gene names for row labels
top_map <- bitr(top_genes, fromType = "ORF",
                toType = "GENENAME", OrgDb = org.Sc.sgd.db)
name_lookup <- setNames(top_map$GENENAME, top_map$ORF)
row_labels <- ifelse(!is.na(name_lookup[rownames(mat)]) & name_lookup[rownames(mat)] != "",
                     name_lookup[rownames(mat)], rownames(mat))
rownames(mat) <- row_labels

annotation_col <- data.frame(
  Stage = sample_table$stage,
  row.names = sample_table$sample_id
)
ann_colors <- list(Stage = c(Early = "#2C5F2D", Thin = "#E8A838", Mature = "#8B2252"))

png(file.path(fig_dir, "figure3_heatmap.png"), width = 8, height = 10, units = "in", res = 300)
pheatmap(mat,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         show_colnames = TRUE,
         fontsize_row = 9,
         fontsize_col = 10,
         main = "Top 30 DE Genes Across Biofilm Stages (LRT)")
dev.off()
cat("Figure 3 (Heatmap) saved.\n")

# ---- Figure 4: GO Biological Process Dotplot ----
if (nrow(as.data.frame(ego_bp)) > 0) {
  go_dotplot <- dotplot(ego_bp, showCategory = 15,
                        title = "GO Biological Process Enrichment (Mature vs Early)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "figure4_go_bp_dotplot.png"), go_dotplot,
         width = 9, height = 7, dpi = 300)
  cat("Figure 4 (GO BP dotplot) saved.\n")
}

# ---- Figure 5: Up vs Down compareCluster ----
if (nrow(as.data.frame(compare_go)) > 0) {
  compare_plot <- dotplot(compare_go, showCategory = 10,
                          title = "GO BP: Upregulated vs Downregulated in Mature Biofilm") +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "figure5_compare_up_down.png"), compare_plot,
         width = 10, height = 8, dpi = 300)
  cat("Figure 5 (Compare clusters) saved.\n")
}

# ---- Figure 6: KEGG Pathway Dotplot ----
if (nrow(as.data.frame(kegg_enrich)) > 0) {
  kegg_dotplot <- dotplot(kegg_enrich, showCategory = 15,
                          title = "KEGG Pathway Enrichment (Mature vs Early)") +
    theme_minimal(base_size = 12)
  ggsave(file.path(fig_dir, "figure6_kegg_dotplot.png"), kegg_dotplot,
         width = 9, height = 7, dpi = 300)
  cat("Figure 6 (KEGG dotplot) saved.\n")
}

# ---- Figure 7: MA Plot ----
png(file.path(fig_dir, "figure7_ma_plot.png"), width = 8, height = 6, units = "in", res = 300)
plotMA(res_early_mature, ylim = c(-5, 5),
       main = "MA Plot: Mature vs Early Biofilm")
dev.off()
cat("Figure 7 (MA plot) saved.\n")

# ---- Summary stats for the report ----
cat("\n=== Summary Statistics ===\n")
cat("Total genes analyzed:", nrow(dds), "\n")
cat("\nEarly vs Thin (padj < 0.05):", sum(res_early_thin$padj < 0.05, na.rm = TRUE), "\n")
cat("Thin vs Mature (padj < 0.05):", sum(res_thin_mature$padj < 0.05, na.rm = TRUE), "\n")
cat("Early vs Mature (padj < 0.05):", sum(res_early_mature$padj < 0.05, na.rm = TRUE), "\n")
cat("LRT (padj < 0.05):", sum(res_lrt$padj < 0.05, na.rm = TRUE), "\n")
cat("\nFigures saved to:", fig_dir, "\n")
