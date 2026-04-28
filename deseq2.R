# -----------------------------
# Load libraries
# -----------------------------
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(org.Mm.eg.db)
library(AnnotationDbi)

# -----------------------------
# Create results directory
# -----------------------------
dir.create("/mnt/d/OJT/AD_multiomics/reference/results", showWarnings = FALSE)

# -----------------------------
# 1. Load featureCounts output
# -----------------------------
counts_raw <- read.table("/mnt/d/OJT/AD_multiomics/reference/raw_counts.txt",
                         header = TRUE,
                         skip = 1,
                         row.names = 1,
                         sep = "\t")

# Extract count columns
counts <- counts_raw[, (ncol(counts_raw)-3):ncol(counts_raw)]
colnames(counts) <- c("ctrl1", "ctrl2", "ad1", "ad2")

# Ensure integer counts
counts <- round(counts)

# -----------------------------
# 2. Metadata
# -----------------------------
coldata <- data.frame(
  row.names = colnames(counts),
  condition = factor(c("control", "control", "AD", "AD"))
)

# Set reference level
coldata$condition <- relevel(coldata$condition, ref = "control")

# -----------------------------
# 3. Create DESeq2 object
# -----------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Filter low counts
dds <- dds[rowSums(counts(dds)) >= 10, ]

# -----------------------------
# 4. Run DESeq2
# -----------------------------
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "AD", "control"))
res <- res[!is.na(res$padj), ]

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

# -----------------------------
# 5. Add gene symbols (NEW)
# -----------------------------
res_df$symbol <- mapIds(org.Mm.eg.db,
                       keys = res_df$gene,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")

# Save full results
write.csv(res_df,
          "/mnt/d/OJT/AD_multiomics/reference/results/DESeq2_all_genes_with_symbols.csv",
          row.names = FALSE)

# -----------------------------
# 6. Significant DEGs
# -----------------------------
sig_degs <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(sig_degs,
          "/mnt/d/OJT/AD_multiomics/reference/results/DESeq2_significant_DEGs.csv",
          row.names = FALSE)

cat("Significant DEGs found:", nrow(sig_degs), "\n")

# -----------------------------
# 7. Volcano Plot
# -----------------------------
res_df$significance <- "Not significant"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Upregulated"
res_df$significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Downregulated"

volcano <- ggplot(res_df,
                  aes(x = log2FoldChange,
                      y = -log10(padj),
                      color = significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "#C0472F",
                                "Downregulated" = "#1B4F8A",
                                "Not significant" = "grey60")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic()

ggsave("/mnt/d/OJT/AD_multiomics/reference/results/volcano_plot.pdf",
       volcano, width = 8, height = 6)

# -----------------------------
# 8. Heatmap (fixed VST issue)
# -----------------------------
if (nrow(dds) < 1000) {
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
} else {
  vsd <- vst(dds, blind = FALSE)
}

topN <- min(50, nrow(sig_degs))
top_genes <- head(sig_degs[order(sig_degs$padj), ], topN)

if (topN > 0) {
  mat <- assay(vsd)[rownames(top_genes), ]
  mat <- mat - rowMeans(mat)

  anno <- data.frame(Condition = coldata$condition)
  rownames(anno) <- colnames(mat)

  pheatmap(mat,
           annotation_col = anno,
           color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100),
           show_rownames = FALSE,
           filename = "/mnt/d/OJT/AD_multiomics/reference/results/heatmap_topDEGs.pdf")
} else {
  cat("No significant genes for heatmap\n")
}
