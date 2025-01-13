## Differential Expression of Genes Analysis


library(DESeq2)
library(pheatmap)
library(ggplot2)

feature_counts <- read.csv("Feature-Counts.csv", row.names = 1, stringsAsFactors = FALSE)

head(feature_counts)

columns_to_remove <- c("Chr", "Start", "End", "Strand", "Length")
feature_counts_clean <- feature_counts[ , !(names(feature_counts) %in% columns_to_remove)]

feature_counts_clean[] <- lapply(feature_counts_clean, as.numeric)

str(feature_counts_clean)

anyNA(feature_counts_clean)

metadata <- data.frame(
  sample = c("SRR12375092", "SRR12375093", "SRR12375094", "SRR12375095", "SRR12375096",
             "SRR12375097", "SRR12375098", "SRR12375099", "SRR12375100", "SRR12375101",
             "SRR12375102", "SRR12375103", "SRR12375104", "SRR12375105", "SRR12375106"),
  GSM = c("GSM4708584", "GSM4708585", "GSM4708586", "GSM4708587", "GSM4708588",
          "GSM4708589", "GSM4708590", "GSM4708591", "GSM4708592", "GSM4708593",
          "GSM4708594", "GSM4708595", "GSM4708596", "GSM4708597", "GSM4708598"),
  condition = c("Control_24h", "Control_24h", "Control_24h", "Control_7d", "Control_7d",
                "Control_7d", "PM_24h_1dose", "PM_24h_1dose", "PM_24h_1dose",
                "PM_24h_30dose", "PM_24h_30dose", "PM_24h_30dose", "PM_7d_1dose",
                "PM_7d_1dose", "PM_7d_1dose"),
  replicate = c(1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3, 1, 2, 3)
)


rownames(metadata) <- metadata$sample

print(metadata)

output_file <- "sample_metadata.csv"
write.csv(metadata, output_file, row.names = FALSE)

cat("Metadata successfully created and saved to:", output_file, "\n")

dim(feature_counts) 
nrow(metadata)       

colnames(feature_counts)

colnames(feature_counts)[!(colnames(feature_counts) %in% rownames(metadata))]

feature_counts <- feature_counts[, rownames(metadata)]

dim(feature_counts)  
nrow(metadata)       

metadata$condition <- factor(metadata$condition)

dds <- DESeqDataSetFromMatrix(countData = feature_counts,
                              colData = metadata,
                              design = ~ condition)

dds

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)


dds <- DESeq(dds)
results <- results(dds)


summary(results)

comparisons <- list(
  "PM_24h_1dose_vs_Control_24h" = c("condition", "PM_24h_1dose", "Control_24h"),
  "PM_24h_30dose_vs_Control_24h" = c("condition", "PM_24h_30dose", "Control_24h"),
  "PM_7d_1dose_vs_Control_7d" = c("condition", "PM_7d_1dose", "Control_7d")
)

deg_summary <- list()

for (name in names(comparisons)) {
  comp <- comparisons[[name]]
  res <- results(dds, contrast = comp)
  
  summary_table <- summary(res)
  upregulated <- sum(res$log2FoldChange > 0 & res$padj < 0.05, na.rm = TRUE)
  downregulated <- sum(res$log2FoldChange < 0 & res$padj < 0.05, na.rm = TRUE)
  
  deg_summary[[name]] <- data.frame(
    Comparison = name,
    Upregulated = upregulated,
    Downregulated = downregulated,
    Total_DEGs = upregulated + downregulated
  )
}

deg_summary_df <- do.call(rbind, deg_summary)
print(deg_summary_df)

deg_summary_df <- data.frame(
  Comparison = c("PM_24h_1dose_vs_Control_24h", "PM_24h_30dose_vs_Control_24h", "PM_7d_1dose_vs_Control_7d"),
  Upregulated = c(18, 1003, 131),
  Downregulated = c(7, 1197, 197)
)

deg_summary_df$Total_DEGs <- deg_summary_df$Upregulated + deg_summary_df$Downregulated


library(reshape2)
deg_melt <- melt(deg_summary_df, id.vars = "Comparison", measure.vars = c("Upregulated", "Downregulated"))

ggplot(deg_melt, aes(x = Comparison, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Upregulated" = "#2b8cbe", "Downregulated" = "#de2d26")) +
  theme_minimal() +
  labs(
    title = "Number of Differentially Expressed Genes",
    x = "Comparison",
    y = "Number of DEGs",
    fill = "Regulation"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

comp_name <- "PM_24h_30dose_vs_Control_24h"
res <- results(dds, contrast = comparisons[[comp_name]])

res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")

ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
  geom_point(alpha = 0.8, size = 1.2) +
  scale_color_manual(values = c("Significant" = "#2b8cbe", "Not Significant" = "#636363")) +
  theme_minimal() +
  labs(
    title = paste("Volcano Plot:", comp_name),
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted P-Value",
    color = "Gene Significance"
  ) +
  theme(legend.position = "top")

top_degs <- rownames(res[order(res$padj, na.last = NA)[1:50], ])

top_counts <- normalized_counts[top_degs, ]

annotation <- data.frame(
  Condition = metadata$condition
)
rownames(annotation) <- rownames(metadata)

library(pheatmap)
pheatmap(
  mat = log2(top_counts + 1), 
  annotation_col = annotation,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  main = "Heatmap of Top 50 DEGs",
  fontsize_row = 6,
  fontsize_col = 8,
  cluster_cols = TRUE,
  cluster_rows = TRUE
)

plotMA(res, ylim = c(-4, 4), main = "MA Plot: PM_24h_30dose vs Control_24h")

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)

ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  theme_minimal() +
  labs(
    title = "PCA of Samples",
    x = paste0("PC1: ", round(100 * attr(pca_data, "percentVar")[1], 1), "% variance"),
    y = paste0("PC2: ", round(100 * attr(pca_data, "percentVar")[2], 1), "% variance"),
    color = "Condition"
  )

generate_volcano <- function(res, comparison_name) {
  res$significant <- ifelse(res$padj < 0.05 & abs(res$log2FoldChange) > 1, "Significant", "Not Significant")
  ggplot(as.data.frame(res), aes(x = log2FoldChange, y = -log10(padj), color = significant)) +
    geom_point(alpha = 0.8, size = 1.2) +
    scale_color_manual(values = c("Significant" = "#1f78b4", "Not Significant" = "#bdbdbd")) +
    theme_minimal() +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value",
      color = "Gene Significance"
    ) +
    theme(legend.position = "top")
}

for (name in names(comparisons)) {
  res <- results(dds, contrast = comparisons[[name]])
  print(generate_volcano(res, name))
}

generate_ma_plot <- function(res, comparison_name) {
  plotMA(res, ylim = c(-4, 4), main = paste("MA Plot:", comparison_name))
}

generate_ma_plot <- function(res, comparison_name) {
  # Check for NA values in adjusted p-value
  if (!any(is.na(res$padj))) {
    plotMA(res, ylim = c(-4, 4), main = paste("MA Plot:", comparison_name))
  } else {
    message(paste("Skipping MA plot for", comparison_name, "due to missing adjusted p-values."))
  }
}

for (name in names(comparisons)) {
  res <- results(dds, contrast = comparisons[[name]])
  generate_ma_plot(res, name)
}

output_file <- "normalized_counts.csv"
write.csv(normalized_counts, file = output_file, row.names = TRUE)

cat("Normalized counts saved to:", output_file, "\n")
