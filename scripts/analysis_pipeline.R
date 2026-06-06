# RNA-seq Differential Expression Analysis of Crohn's Disease
# Author: Hazrat Maghaz
# Workflow: DESeq2 differential expression analysis + visualization + Hallmark GSEA

# -----------------------------
# 1. Load required packages
# -----------------------------

required_packages <- c(
  "DESeq2",
  "tidyverse",
  "dplyr",
  "ggplot2",
  "pheatmap",
  "EnhancedVolcano",
  "RColorBrewer",
  "clusterProfiler",
  "msigdbr",
  "magrittr"
)

missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  stop(
    "The following packages are missing: ",
    paste(missing_packages, collapse = ", "),
    "\nPlease install them before running this script."
  )
}

library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(RColorBrewer)
library(clusterProfiler)
library(msigdbr)
library(magrittr)

# -----------------------------
# 2. Define project paths
# -----------------------------

data_dir <- "data"
figures_dir <- file.path("results", "figures")
tables_dir <- file.path("results", "tables")

dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

counts_file <- file.path(data_dir, "Counts.tsv")
metadata_file <- file.path(data_dir, "Metadata.tsv")

# -----------------------------
# 3. Load input data
# -----------------------------

metadata <- read.delim(metadata_file, check.names = FALSE)
counts <- read.delim(counts_file, row.names = 1, check.names = FALSE)

# Use the first metadata column as sample IDs if row names are not already set
if (!all(colnames(counts) %in% rownames(metadata))) {
  rownames(metadata) <- metadata[[1]]
}

metadata <- metadata[colnames(counts), , drop = FALSE]

if (!"Category" %in% colnames(metadata)) {
  stop("Metadata must contain a 'Category' column for disease/control grouping.")
}

metadata$Category <- as.factor(metadata$Category)

# -----------------------------
# 4. Create DESeq2 dataset
# -----------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ Category
)

# Remove very low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Set control group as reference if present
control_label <- "non inflammatory bowel disease control"

if (control_label %in% levels(dds$Category)) {
  dds$Category <- relevel(dds$Category, ref = control_label)
}

# -----------------------------
# 5. Run DESeq2 analysis
# -----------------------------

dds <- DESeq(dds)
res <- results(dds, alpha = 0.05)
res_ordered <- res[order(res$padj), ]

write.csv(
  as.data.frame(res_ordered),
  file = file.path(tables_dir, "DESeq2_results_ordered.csv")
)

res_sig <- subset(as.data.frame(res_ordered), padj < 0.05)

write.csv(
  res_sig,
  file = file.path(tables_dir, "DESeq2_significant_DEGs_FDR_0.05.csv")
)

# -----------------------------
# 6. Variance-stabilizing transformation
# -----------------------------

rld <- vst(dds, blind = FALSE)

# -----------------------------
# 7. Sample distance heatmap
# -----------------------------

sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)

colors <- colorRampPalette(rev(brewer.pal(9, "Greens")))(255)

png(file.path(figures_dir, "Heatmap.png"), width = 1200, height = 1000, res = 150)
pheatmap(
  sample_dist_matrix,
  clustering_distance_rows = sample_dists,
  clustering_distance_cols = sample_dists,
  col = colors,
  main = "Sample-to-sample distance heatmap"
)
dev.off()

# -----------------------------
# 8. PCA plot
# -----------------------------

pca_plot <- plotPCA(rld, intgroup = "Category") +
  ggtitle("PCA of Crohn's Disease and Control Samples") +
  theme_minimal()

ggsave(
  filename = file.path(figures_dir, "PCA_plot.png"),
  plot = pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

# -----------------------------
# 9. Volcano plot
# -----------------------------

png(file.path(figures_dir, "Volcano_plot.png"), width = 1200, height = 1000, res = 150)
EnhancedVolcano(
  res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "padj",
  title = "Crohn's Disease vs Control",
  subtitle = "DESeq2 differential expression analysis",
  pCutoff = 0.05,
  FCcutoff = 1.0,
  pointSize = 2.0,
  labSize = 3.0
)
dev.off()

# -----------------------------
# 10. Hallmark GSEA
# -----------------------------

res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)

ranked_genes <- res_df$stat
names(ranked_genes) <- res_df$gene
ranked_genes <- ranked_genes[!is.na(ranked_genes)]
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

hallmark_sets <- msigdbr(species = "Homo sapiens", category = "H") %>%
  dplyr::select(gs_name, gene_symbol)

gsea_results <- GSEA(
  geneList = ranked_genes,
  TERM2GENE = hallmark_sets,
  pvalueCutoff = 0.05,
  verbose = FALSE
)

gsea_table <- as.data.frame(gsea_results)

write.csv(
  gsea_table,
  file = file.path(tables_dir, "Hallmark_GSEA_results.csv"),
  row.names = FALSE
)

if (nrow(gsea_table) > 0) {
  png(file.path(figures_dir, "Hallmark_plot.png"), width = 1200, height = 900, res = 150)
  print(dotplot(gsea_results, showCategory = 15) + ggtitle("Hallmark GSEA Results"))
  dev.off()
}

# -----------------------------
# 11. Save session information
# -----------------------------

sink(file.path(tables_dir, "sessionInfo.txt"))
sessionInfo()
sink()

message("Analysis completed successfully.")
message("Figures saved in: ", figures_dir)
message("Tables saved in: ", tables_dir)
