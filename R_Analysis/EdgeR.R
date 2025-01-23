#############################Edge_R_Salmon_on_Raw_Data######################

library(edgeR)
library(ggplot2)


path_to_counts <- "Salmon_merged_quant_counts.csv"
path_to_metadata <- "metadata.csv"


counts <- read.csv(path_to_counts, row.names = 1)
metadata <- read.csv(path_to_metadata)


counts <- round(as.matrix(counts))

metadata <- metadata[match(colnames(counts), metadata$Run), ]
row.names(metadata) <- metadata$Run
metadata$Run <- NULL

metadata$treatment <- as.factor(metadata$treatment)
levels(metadata$treatment) <- make.names(levels(metadata$treatment))  # Ensure valid names


dge <- DGEList(counts = counts, group = metadata$treatment)


keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~ 0 + metadata$treatment)
colnames(design) <- levels(metadata$treatment)

dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)

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


create_volcano_plot <- function(res_table, condition1, condition2) {
  res_table$Significant <- res_table$FDR < 0.05 & abs(res_table$logFC) > 1
  ggplot(res_table, aes(x = logFC, y = -log10(PValue), color = Significant)) +
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
}


for (comparison in comparisons) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  

  contrast <- makeContrasts(
    contrasts = paste0(condition2, "-", condition1),
    levels = design
  )
  

  qlf <- glmQLFTest(fit, contrast = contrast)
  res_table <- topTags(qlf, n = Inf)$table  # Get all genes

  upregulated <- sum(res_table$logFC > 1 & res_table$FDR < 0.05, na.rm = TRUE)
  downregulated <- sum(res_table$logFC < -1 & res_table$FDR < 0.05, na.rm = TRUE)
  

  bar_data <- rbind(
    bar_data,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
  
  volcano_plot <- create_volcano_plot(res_table, condition1, condition2)
  print(volcano_plot)  # Display plot in RStudio
  ggsave(paste0("Volcano_", condition1, "_vs_", condition2, ".png"), volcano_plot, width = 8, height = 6)
}


bar_plot <- ggplot(bar_data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), vjust = -0.5, 
            position = position_dodge(width = 0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Upregulated and Downregulated Genes per Comparison (edgeR)",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#1f78b4", "Downregulated" = "#e31a1c"))

print(bar_plot)
ggsave("DE_Gene_Barplot.png", bar_plot, width = 10, height = 6)

#########################EdgeR_Salmon_on_Filtered_Data_######################################

path_to_filtered_counts <- "Salmon_Filtered_merged_quant_counts.csv"  
path_to_metadata <- "metadata.csv"

filtered_counts <- read.csv(path_to_filtered_counts, row.names = 1) 
metadata <- read.csv(path_to_metadata)

filtered_counts <- round(as.matrix(filtered_counts))  


metadata <- metadata[match(colnames(filtered_counts), metadata$Run), ] 
row.names(metadata) <- metadata$Run
metadata$Run <- NULL


metadata$treatment <- as.factor(metadata$treatment)
levels(metadata$treatment) <- make.names(levels(metadata$treatment))

dge <- DGEList(counts = filtered_counts, group = metadata$treatment)  

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]


dge <- calcNormFactors(dge)


design <- model.matrix(~0 + metadata$treatment)
colnames(design) <- levels(metadata$treatment)


dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)


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


create_volcano_plot <- function(res_table, condition1, condition2) {
  res_table$Significant <- res_table$FDR < 0.05 & abs(res_table$logFC) > 1
  ggplot(res_table, aes(x = logFC, y = -log10(PValue), color = Significant)) +
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
}

for (comparison in comparisons) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  contrast <- makeContrasts(
    contrasts = paste0(condition2, "-", condition1),
    levels = design
  )
  
  qlf <- glmQLFTest(fit, contrast = contrast)
  res_table <- topTags(qlf, n = Inf)$table
  
  upregulated <- sum(res_table$logFC > 1 & res_table$FDR < 0.05, na.rm = TRUE)
  downregulated <- sum(res_table$logFC < -1 & res_table$FDR < 0.05, na.rm = TRUE)
  
  bar_data <- rbind(
    bar_data,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
  
  volcano_plot <- create_volcano_plot(res_table, condition1, condition2)
  print(volcano_plot)
  ggsave(paste0("Volcano_", condition1, "_vs_", condition2, ".png"), volcano_plot, width = 8, height = 6)
}


bar_plot <- ggplot(bar_data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), vjust = -0.5, 
            position = position_dodge(width = 0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Upregulated and Downregulated Genes per Comparison (edgeR)",
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#1f78b4", "Downregulated" = "#e31a1c"))


print(bar_plot)
#ggsave("DE_Gene_Barplot_Filtered.png", bar_plot, width = 10, height = 6)  # Updated output name

#####################################EdgeR_Kallisto_on_Raw_Data##################################

path_to_kallisto_counts <- "Kallisto_merged_kallisto_counts.csv" 
path_to_metadata <- "metadata.csv"

kallisto_counts <- read.csv(path_to_kallisto_counts, row.names = 1)  
metadata <- read.csv(path_to_metadata)

kallisto_counts <- round(as.matrix(kallisto_counts)) 

metadata <- metadata[match(colnames(kallisto_counts), metadata$Run), ]  
row.names(metadata) <- metadata$Run
metadata$Run <- NULL

metadata$treatment <- as.factor(metadata$treatment)
levels(metadata$treatment) <- make.names(levels(metadata$treatment))

dge <- DGEList(counts = kallisto_counts, group = metadata$treatment)  


keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~0 + metadata$treatment)
colnames(design) <- levels(metadata$treatment)

dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)

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

create_volcano_plot <- function(res_table, condition1, condition2) {
  res_table$Significant <- res_table$FDR < 0.05 & abs(res_table$logFC) > 1
  ggplot(res_table, aes(x = logFC, y = -log10(PValue), color = Significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    theme_minimal() +
    scale_color_manual(values = c("gray", "blue")) +  # Magenta for significant
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
}


for (comparison in comparisons) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  contrast <- makeContrasts(
    contrasts = paste0(condition2, "-", condition1),
    levels = design
  )
  
  qlf <- glmQLFTest(fit, contrast = contrast)
  res_table <- topTags(qlf, n = Inf)$table
  
  upregulated <- sum(res_table$logFC > 1 & res_table$FDR < 0.05, na.rm = TRUE)
  downregulated <- sum(res_table$logFC < -1 & res_table$FDR < 0.05, na.rm = TRUE)
  
  bar_data <- rbind(
    bar_data,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
  
  volcano_plot <- create_volcano_plot(res_table, condition1, condition2)
  print(volcano_plot)
  ggsave(paste0("Volcano_Kallisto_", condition1, "_vs_", condition2, ".png"), volcano_plot, width = 8, height = 6)  # Updated prefix
}

# Bar plot with new colors (purple/orange)
bar_plot <- ggplot(bar_data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), vjust = -0.5, 
            position = position_dodge(width = 0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Up/Downregulated Genes: Kallisto Dataset",  # Updated title
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#87CEEB", "Downregulated" = "#FFA07A"))  # Purple/Orange

print(bar_plot)
ggsave("DE_Gene_Barplot_Kallisto.png", bar_plot, width = 10, height = 6)  # Updated output name

############################ EdgeR_Kallisto_Filtered_Data ########################

path_to_kallisto_filtered_counts <- "Kallisto_filtered_merged_kallisto_counts.csv"  
path_to_metadata <- "metadata.csv"

kallisto_filtered_counts <- read.csv(path_to_kallisto_filtered_counts, row.names = 1)  
metadata <- read.csv(path_to_metadata)

kallisto_filtered_counts <- round(as.matrix(kallisto_filtered_counts))  


metadata <- metadata[match(colnames(kallisto_filtered_counts), metadata$Run), ]  
row.names(metadata) <- metadata$Run
metadata$Run <- NULL

metadata$treatment <- as.factor(metadata$treatment)
levels(metadata$treatment) <- make.names(levels(metadata$treatment))

dge <- DGEList(counts = kallisto_filtered_counts, group = metadata$treatment)  

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

design <- model.matrix(~0 + metadata$treatment)
colnames(design) <- levels(metadata$treatment)

dge <- estimateDisp(dge, design)

fit <- glmQLFit(dge, design)

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

create_volcano_plot <- function(res_table, condition1, condition2) {
  res_table$Significant <- res_table$FDR < 0.05 & abs(res_table$logFC) > 1
  ggplot(res_table, aes(x = logFC, y = -log10(PValue), color = Significant)) +
    geom_point(alpha = 0.5, size = 1.2) +
    theme_minimal() +
    scale_color_manual(values = c("gray", "blue")) +  # Dark blue for significant
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
}

for (comparison in comparisons) {
  condition1 <- comparison[1]
  condition2 <- comparison[2]
  
  contrast <- makeContrasts(
    contrasts = paste0(condition2, "-", condition1),
    levels = design
  )
  
  qlf <- glmQLFTest(fit, contrast = contrast)
  res_table <- topTags(qlf, n = Inf)$table
  
  upregulated <- sum(res_table$logFC > 1 & res_table$FDR < 0.05, na.rm = TRUE)
  downregulated <- sum(res_table$logFC < -1 & res_table$FDR < 0.05, na.rm = TRUE)
  
  bar_data <- rbind(
    bar_data,
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Upregulated",
               Count = upregulated),
    data.frame(Comparison = paste(condition1, "vs", condition2),
               Direction = "Downregulated",
               Count = downregulated)
  )
  
  volcano_plot <- create_volcano_plot(res_table, condition1, condition2)
  print(volcano_plot)
  ggsave(paste0("Volcano_Kallisto_Filtered_", condition1, "_vs_", condition2, ".png"), 
         volcano_plot, width = 8, height = 6)  # Updated prefix
}

bar_plot <- ggplot(bar_data, aes(x = Comparison, y = Count, fill = Direction)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = Count), vjust = -0.5, 
            position = position_dodge(width = 0.9), size = 4) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold")
  ) +
  labs(
    title = "Up/Downregulated Genes: Kallisto Filtered Dataset",  # Updated title
    x = "Comparison",
    y = "Number of Genes",
    fill = "Direction"
  ) +
  scale_fill_manual(values = c("Upregulated" = "#87CEEB", "Downregulated" = "#FFA07A"))  # Dark green/gold

print(bar_plot)
ggsave("DE_Gene_Barplot_Kallisto_Filtered.png", bar_plot, width = 10, height = 6)  # Updated output name