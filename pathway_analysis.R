# ── pathway_analysis.R (gprofiler2 version) ─────────────────────────────────
# Replaces clusterProfiler/enrichplot with gprofiler2
# No C++ compilation required — queries Ensembl API online
# Compatible with Ensembl IDs directly (no bitr() conversion needed)

library(gprofiler2)
library(ggplot2)

# ── Install if needed ────────────────────────────────────────────────────────
# install.packages("gprofiler2")

# ── Load DEGs ────────────────────────────────────────────────────────────────
sig_degs <- read.csv("results/DESeq2_significant_DEGs.csv")

cat("Loaded", nrow(sig_degs), "significant DEGs\n")
cat("Gene ID preview:", head(sig_degs$gene, 3), "\n")

# ── Run GO + KEGG enrichment in one call ─────────────────────────────────────
results <- gost(
  query             = sig_degs$gene,
  organism          = "mmusculus",        # Mus musculus
  sources           = c("GO:BP", "GO:MF", "GO:CC", "KEGG", "REAC"),
  correction_method = "fdr",
  significant       = TRUE,
  evcodes           = TRUE               # include evidence codes
)

if (is.null(results)) {
  stop("No significant pathways found. Check gene IDs or relax FDR threshold.")
}

cat("Enrichment complete.", nrow(results$result), "significant terms found\n")

# ── Interactive Manhattan plot (opens in browser) ─────────────────────────────
p_interactive <- gostplot(results, interactive = TRUE)
# Uncomment to view: p_interactive

# ── Static publication plot ───────────────────────────────────────────────────
p_static <- gostplot(results, interactive = FALSE)
ggsave("results/gprofiler_manhattan.pdf", plot = p_static, width = 12, height = 7)
cat("Saved: results/gprofiler_manhattan.pdf\n")

# ── Split results by source ───────────────────────────────────────────────────
res_df <- results$result

go_bp   <- res_df[res_df$source == "GO:BP", ]
go_mf   <- res_df[res_df$source == "GO:MF", ]
go_cc   <- res_df[res_df$source == "GO:CC", ]
kegg    <- res_df[res_df$source == "KEGG",  ]
reactome <- res_df[res_df$source == "REAC", ]

cat("\nTerms per source:\n")
cat("  GO:BP    :", nrow(go_bp), "\n")
cat("  GO:MF    :", nrow(go_mf), "\n")
cat("  GO:CC    :", nrow(go_cc), "\n")
cat("  KEGG     :", nrow(kegg),  "\n")
cat("  Reactome :", nrow(reactome), "\n")

# ── Dotplot for GO:BP (top 20) ────────────────────────────────────────────────
if (nrow(go_bp) > 0) {
  go_bp_top <- head(go_bp[order(go_bp$p_value), ], 20)

  p_gobp <- ggplot(go_bp_top, aes(
    x    = -log10(p_value),
    y    = reorder(term_name, -p_value),
    size = intersection_size,
    color = p_value
  )) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    scale_size_continuous(range = c(3, 10)) +
    labs(
      title  = "GO Biological Process - Top 20",
      x      = "-log10(FDR)",
      y      = NULL,
      size   = "Gene count",
      color  = "FDR p-value"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))

  ggsave("results/GO_BP_dotplot.pdf", plot = p_gobp, width = 10, height = 8)
  cat("Saved: results/GO_BP_dotplot.pdf\n")
}

# ── Dotplot for KEGG (top 15) ─────────────────────────────────────────────────
if (nrow(kegg) > 0) {
  kegg_top <- head(kegg[order(kegg$p_value), ], 15)

  p_kegg <- ggplot(kegg_top, aes(
    x     = -log10(p_value),
    y     = reorder(term_name, -p_value),
    size  = intersection_size,
    color = p_value
  )) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    scale_size_continuous(range = c(3, 10)) +
    labs(
      title = "KEGG Pathway Enrichment — Top 15",
      x     = "-log10(FDR)",
      y     = NULL,
      size  = "Gene count",
      color = "FDR p-value"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.y = element_text(size = 9))

  ggsave("results/KEGG_dotplot.pdf", plot = p_kegg, width = 10, height = 7)
  cat("Saved: results/KEGG_dotplot.pdf\n")
}

# ── Flatten list columns before saving ───────────────────────────────────────
flatten_df <- function(df) {
  for (col in names(df)) {
    if (is.list(df[[col]])) {
      df[[col]] <- sapply(df[[col]], function(x) paste(x, collapse = ","))
    }
  }
  return(df)
}

res_flat      <- flatten_df(res_df)
go_bp_flat    <- flatten_df(go_bp)
go_mf_flat    <- flatten_df(go_mf)
go_cc_flat    <- flatten_df(go_cc)
kegg_flat     <- flatten_df(kegg)
reactome_flat <- flatten_df(reactome)

write.csv(go_bp_flat,    "results/GO_BP_results.csv",        row.names = FALSE)
write.csv(go_mf_flat,    "results/GO_MF_results.csv",        row.names = FALSE)
write.csv(go_cc_flat,    "results/GO_CC_results.csv",        row.names = FALSE)
write.csv(kegg_flat,     "results/KEGG_results.csv",         row.names = FALSE)
write.csv(reactome_flat, "results/Reactome_results.csv",     row.names = FALSE)
write.csv(res_flat,      "results/gprofiler_all_results.csv", row.names = FALSE)

cat("\nAll results saved to results/\n")
cat("Pathway analysis complete.\n")
cat("Carry forward results/DESeq2_significant_DEGs.csv to ChIP-seq/ChIPseeker step.\n")
