# Load necessary libraries
library(DESeq2)
library(readr)
library(bigPint)

path_to_counts <- "Salmon_merged_quant_counts.csv"
path_to_metadata <- "metadata.csv"

counts <- read.csv(path_to_counts, row.names = 1)
metadata <- read.csv(path_to_metadata)

counts <- round(as.matrix(counts))

if (!is.matrix(counts)) stop("Counts data is not a matrix.")
if (!all(counts == round(counts))) stop("Counts data contains non-integer values.")

if (!"Run" %in% colnames(metadata)) stop("'Run' column is missing in metadata.")
if (!"treatment" %in% colnames(metadata)) stop("'treatment' column is missing in metadata.")

metadata <- metadata[match(colnames(counts), metadata$Run), ]

if (!all(colnames(counts) == metadata$Run)) {
  stop("Metadata and count data columns are not aligned. Check the metadata or counts file.")
} else {
  cat("Metadata and count data are properly aligned.\n")
}

row.names(metadata) <- metadata$Run
metadata$Run <- NULL

metadata$treatment <- as.factor(metadata$treatment)

dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData = metadata,
  design = ~ treatment
)

keep <- rowSums(counts(dds)) > 10
dds <- dds[keep, ]

dds <- DESeq(dds)

res <- results(dds)

levels(dds$treatment)

comparisons <- list(
  c("PM_24h_1dose", "PM_24h_30dose"),
  c("Control_24h", "PM_24h_30dose"),
  c("Control_7d", "PM_7d_1dose"),
  c("Control_24h", "PM_24h_1dose")
)


create_volcano_plot <- function(res, condition1, condition2) {
  res_df <- as.data.frame(res)
  res_df$Significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    theme_minimal() +
    scale_color_manual(values = c("gray", "red")) +
    labs(
      title = paste("Volcano Plot:", condition1, "vs", condition2),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "top"
    )
  
  print(p)
}

for (comparison in comparisons) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]

  res <- results(dds, contrast = c("treatment", condition2, condition1))

  create_volcano_plot(res, condition1, condition2)
}


comparisons <- list(
  c("PM_24h_1dose", "PM_24h_30dose"),
  c("Control_24h", "PM_24h_30dose"),
  c("Control_7d", "PM_7d_1dose"),
  c("Control_24h", "PM_24h_1dose")
)

bar_data <- data.frame(
  Comparison = character(),
  Direction = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

for (comparison in comparisons) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]

  res <- results(dds, contrast = c("treatment", condition2, condition1))

  upregulated <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
  downregulated <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)

  bar_data <- rbind(
    bar_data,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
}



ggplot(bar_data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Upregulated and Downregulated Genes per Comparison",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#1f78b4", "Downregulated" = "#e31a1c"))

################################################################ DESeq2_filterd_Salmon ###############################################################################

path_to_filtered_counts <- "Salmon_Filtered_merged_quant_counts.csv"
path_to_metadata <- "metadata.csv"

filtered_counts <- read.csv(path_to_filtered_counts, row.names = 1)
metadata <- read.csv(path_to_metadata)

filtered_counts <- round(as.matrix(filtered_counts))

metadata <- metadata[match(colnames(filtered_counts), metadata$Run), ]
row.names(metadata) <- metadata$Run
metadata$Run <- NULL
metadata$treatment <- as.factor(metadata$treatment)

dds_filtered <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = metadata,
  design = ~ treatment
)

keep_filtered <- rowSums(counts(dds_filtered)) > 10
dds_filtered <- dds_filtered[keep_filtered, ]

dds_filtered <- DESeq(dds_filtered)

comparisons_filtered <- list(
  c("PM_24h_1dose", "PM_24h_30dose"),
  c("Control_24h", "PM_24h_30dose"),
  c("Control_7d", "PM_7d_1dose"),
  c("Control_24h", "PM_24h_1dose")
)

create_volcano_plot <- function(res, condition1, condition2) {
  res_df <- as.data.frame(res)
  res_df$Significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  

  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    theme_minimal() +
    scale_color_manual(values = c("gray", "red")) +
    labs(
      title = paste("Volcano Plot:", condition1, "vs", condition2),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "top"
    )
  
  print(p)
}

for (comparison in comparisons_filtered) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  res <- results(dds_filtered, contrast = c("treatment", condition2, condition1))
  
  create_volcano_plot(res, condition1, condition2)
}

bar_data_filtered <- data.frame(
  Comparison = character(),
  Direction = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

for (comparison in comparisons_filtered) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  res <- results(dds_filtered, contrast = c("treatment", condition2, condition1))
  
  upregulated <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
  downregulated <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)
  
  bar_data_filtered <- rbind(
    bar_data_filtered,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
}

ggplot(bar_data_filtered, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Upregulated and Downregulated Genes per Comparison (Filtered Data)",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#1f78b4", "Downregulated" = "#e31a1c"))


###################################################### DESeq2_Kallisto_Count_File############################################################

path_to_kallisto_counts <- "Kallisto_merged_kallisto_counts.csv"
path_to_metadata <- "metadata.csv"

kallisto_counts <- read.csv(path_to_kallisto_counts, row.names = 1)
metadata <- read.csv(path_to_metadata)

kallisto_counts <- round(as.matrix(kallisto_counts))

metadata <- metadata[match(colnames(kallisto_counts), metadata$Run), ]
row.names(metadata) <- metadata$Run
metadata$Run <- NULL
metadata$treatment <- as.factor(metadata$treatment)

dds_kallisto <- DESeqDataSetFromMatrix(
  countData = kallisto_counts,
  colData = metadata,
  design = ~ treatment
)

keep_kallisto <- rowSums(counts(dds_kallisto)) > 10
dds_kallisto <- dds_kallisto[keep_kallisto, ]

dds_kallisto <- DESeq(dds_kallisto)

comparisons_kallisto <- list(
  c("PM_24h_1dose", "PM_24h_30dose"),
  c("Control_24h", "PM_24h_30dose"),
  c("Control_7d", "PM_7d_1dose"),
  c("Control_24h", "PM_24h_1dose")
)

create_volcano_plot <- function(res, condition1, condition2) {
  res_df <- as.data.frame(res)
  res_df$Significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    theme_minimal() +
    scale_color_manual(values = c("gray", "blue")) +
    labs(
      title = paste("Volcano Plot (Kallisto):", condition1, "vs", condition2),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "top"
    )
  
  print(p)
}

for (comparison in comparisons_kallisto) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  res <- results(dds_kallisto, contrast = c("treatment", condition2, condition1))
  
  create_volcano_plot(res, condition1, condition2)
}

bar_data_kallisto <- data.frame(
  Comparison = character(),
  Direction = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

for (comparison in comparisons_kallisto) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  res <- results(dds_kallisto, contrast = c("treatment", condition2, condition1))
  
  upregulated <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
  downregulated <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)
  
  bar_data_kallisto <- rbind(
    bar_data_kallisto,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
}

ggplot(bar_data_kallisto, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Upregulated and Downregulated Genes per Comparison (Kallisto Data)",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#87CEEB", "Downregulated" = "#FFA07A"))

####################################################################### DESeq2_Kallisto_Filtered_##########################################################################

path_to_kallisto_filtered_counts <- "Kallisto_filtered_merged_kallisto_counts.csv"
path_to_metadata <- "metadata.csv"

kallisto_filtered_counts <- read.csv(path_to_kallisto_filtered_counts, row.names = 1)
metadata <- read.csv(path_to_metadata)

kallisto_filtered_counts <- round(as.matrix(kallisto_filtered_counts))

metadata <- metadata[match(colnames(kallisto_filtered_counts), metadata$Run), ]
row.names(metadata) <- metadata$Run
metadata$Run <- NULL
metadata$treatment <- as.factor(metadata$treatment)

dds_kallisto_filtered <- DESeqDataSetFromMatrix(
  countData = kallisto_filtered_counts,
  colData = metadata,
  design = ~ treatment
)

keep_kallisto_filtered <- rowSums(counts(dds_kallisto_filtered)) > 10
dds_kallisto_filtered <- dds_kallisto_filtered[keep_kallisto_filtered, ]

dds_kallisto_filtered <- DESeq(dds_kallisto_filtered)

comparisons_kallisto_filtered <- list(
  c("PM_24h_1dose", "PM_24h_30dose"),
  c("Control_24h", "PM_24h_30dose"),
  c("Control_7d", "PM_7d_1dose"),
  c("Control_24h", "PM_24h_1dose")
)

create_volcano_plot <- function(res, condition1, condition2) {
  res_df <- as.data.frame(res)
  res_df$Significant <- res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1
  
  p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = Significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    theme_minimal() +
    scale_color_manual(values = c("gray", "blue")) +
    labs(
      title = paste("Volcano Plot (Kallisto Filtered):", condition1, "vs", condition2),
      x = "Log2 Fold Change",
      y = "-Log10 P-value"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = "top"
    )
  
  print(p)
}

for (comparison in comparisons_kallisto_filtered) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  res <- results(dds_kallisto_filtered, contrast = c("treatment", condition2, condition1))
  
  create_volcano_plot(res, condition1, condition2)
}

bar_data_kallisto_filtered <- data.frame(
  Comparison = character(),
  Direction = character(),
  Count = integer(),
  stringsAsFactors = FALSE
)

for (comparison in comparisons_kallisto_filtered) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  res <- results(dds_kallisto_filtered, contrast = c("treatment", condition2, condition1))
  
  upregulated <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
  downregulated <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)
  
  bar_data_kallisto_filtered <- rbind(
    bar_data_kallisto_filtered,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
}

ggplot(bar_data_kallisto_filtered, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_text(aes(label = Count), vjust = -0.5, position = position_dodge(0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Upregulated and Downregulated Genes per Comparison (Kallisto Filtered Data)",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#87CEEB", "Downregulated" = "#FFA07A"))
