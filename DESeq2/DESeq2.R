# Step 1: Load the necessary Libraries and the Gene_Counts_Data
#--------------------------------------------------

# Load dplyr for data manipulation
library(dplyr)

# Load DESeq2 for differential gene expression analysis
library(DESeq2) 

# Load clusterProfiler and org.Hs.eg.db for enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)

# Load ggplot2 for plotting and EnhancedVolcano for visualizing DEGs
library(ggplot2)
library(EnhancedVolcano)


# Gene counts data. Ensure you set the correct working directory or provide the full path to the file.
counts <- read.csv("gene_counts_only.csv", row.names = 1)

# View first few lines
head(counts)

# Step 2: Create Sample Information
#-----------------------------------
# Define metadata for each sample, including experimental conditions and time points.
# Ensure the sample identifiers in row.names match the column names in 'counts' data.

sample_info <- data.frame(
  row.names = c("SRR22269883", "SRR22269882", "SRR22269881", 
                "SRR22269880", "SRR22269879", "SRR22269878", 
                "SRR22269877", "SRR22269876", "SRR22269875", 
                "SRR22269874", "SRR22269873", "SRR22269872"),
  Condition = c(rep("Mock", 3), rep("SARS_CoV_2", 3), rep("Mock", 3), rep("SARS_CoV_2", 3)),
  TimePoint = c(rep(24, 6), rep(72, 6))
)
sample_info$Condition <- factor(gsub("-", "_", sample_info$Condition))
sample_info$TimePoint <- factor(paste0("T", sample_info$TimePoint))

# Match the order of columns in counts with sample_info rows
counts <- counts[, rownames(sample_info)]

# Verify alignment between counts data and sample_info
head(sample_info)

# Step 3: DESeq2 Analysis
#-------------------------

# Create DESeq2 dataset object from count matrix and sample metadata
# Specify design formula to include both Condition and TimePoint factors
dds <- DESeqDataSetFromMatrix(countData = counts, colData = sample_info, design = ~ Condition + TimePoint)
dds <- DESeq(dds)

# Differential expression analysis: 72H vs 24H
# --------------------------------------------
# Extract results comparing TimePoint T72 to T24
res_time <- results(dds, contrast = c("TimePoint", "T72", "T24"))

# Remove any rows with NA values in adjusted p-values to avoid errors in further filtering
res_time <- res_time[!is.na(res_time$padj), ]

# Filter significant DEGs (differentially expressed genes) based on adjusted p-value threshold
degs_time <- res_time[res_time$padj < 0.2, ]

# Display summary of results for the time-based comparison
summary(degs_time)

# Differential expression analysis:SARS-CoV-2 vs Mock
#----------------------------------------------------
# Extract results comparing Condition SARS_CoV_2 to Mock
res_condition <- results(dds, contrast = c("Condition", "SARS_CoV_2", "Mock"))

# Remove any rows with NA values in adjusted p-values
res_condition <- res_condition[!is.na(res_condition$padj), ]


# Filter significant DEGs for condition-based comparison
degs_condition <- res_condition[res_condition$padj < 0.2, ]

# Display summary of results for the condition-based comparison
summary(degs_condition)

# Step 4: Visualization of the DESQ2 results
#----------------------------------------------------------------

# MA Plot for TimePoint comparison
plotMA(res_time, ylim = c(-5, 5), main = "MA Plot: 72H vs 24H", cex = 1.0)

# MA Plot for Condition comparison
plotMA(res_condition, ylim = c(-5, 5), main = "MA Plot: SARS-CoV-2 vs Mock", cex = 1.0)

#---------------------------------------------------------------

# Volcano Plot for TimePoint comparison
EnhancedVolcano(res_time,
                lab = rownames(res_time),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot: 72H vs 24H',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3.0)

# Volcano Plot for SARS-CoV-2 vs Mock comparison.
EnhancedVolcano(res_condition,
                lab = rownames(res_condition),
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Volcano Plot: SARS-CoV-2 vs Mock',
                pCutoff = 0.05,                      
                FCcutoff = 1,                        
                pointSize = 3.0,                     
                labSize = 4.0,                       
                xlim = c(-3, 3),                     
                ylim = c(0, -log10(min(res_condition$pvalue, na.rm = TRUE)) + 2),
                col = c('grey', 'skyblue', 'orange', 'red'),  
                legendPosition = "top",
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,               
                widthConnectors = 0.5,
                colConnectors = 'black',
                gridlines.major = FALSE,
                gridlines.minor = FALSE)

#-----------------------------------------------------------------
# Variance stabilizing transformation for PCA
vsd <- vst(dds, blind = FALSE)  

# Extract PCA data
pca_data <- plotPCA(vsd, intgroup = c("Condition", "TimePoint"), returnData = TRUE)

# Define a color palette and shapes for each condition
color_palette <- c("Mock" = "#1f77b4", "SARS_CoV_2" = "#ff7f0e")
shape_values <- c("T24" = 16, "T72" = 17)  # Different shapes for time points (circle and triangle)


ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = TimePoint)) +
  geom_point(size = 4, alpha = 0.8) +                          
  scale_color_manual(values = color_palette) +                 
  scale_shape_manual(values = shape_values) +                  
  labs(title = "PCA Plot: Condition and TimePoint", 
       x = paste0("PC1: ", round(attr(pca_data, "percentVar")[1] * 100), "% Variance"),
       y = paste0("PC2: ", round(attr(pca_data, "percentVar")[2] * 100), "% Variance")) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 14),
    legend.position = "right",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )
#-----------------------------------------------------------------

# Step 5: GO Enrichment Analysis
#-------------------------------
# Extract list of significant genes for GO enrichment analysis.
degs_time_genes <- rownames(degs_time)
degs_condition_genes <- rownames(degs_condition)

## Perform GO enrichment analysis for time-based comparison (72H vs 24H)
go_enrichment_time <- enrichGO(gene = degs_time_genes, OrgDb = org.Hs.eg.db,
                               keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                               pvalueCutoff = 0.1, qvalueCutoff = 0.2,
                               minGSSize = 5, maxGSSize = 500, readable = TRUE)

head(go_enrichment_time_df)

# Perform GO enrichment analysis for condition-based comparison (SARS-CoV-2 vs Mock)
go_enrichment_condition <- enrichGO(gene = degs_condition_genes, OrgDb = org.Hs.eg.db,
                                    keyType = "ENSEMBL", ont = "BP", pAdjustMethod = "BH",
                                    pvalueCutoff = 0.5, qvalueCutoff = 0.4,
                                    minGSSize = 2, maxGSSize = 2000, readable = TRUE)

head(go_enrichment_condition_df)

# Step 6: Visualization of GO Enrichment
#---------------------------------------
# Convert enrichment results to data frames
go_enrichment_time_df <- as.data.frame(go_enrichment_time)
go_enrichment_condition_df <- as.data.frame(go_enrichment_condition)

# Select the top 20 GO terms by Fold Enrichment for time-based comparison
go_enrichment_time_top <- go_enrichment_time_df %>%
  arrange(desc(FoldEnrichment)) %>%
  head(20)

# Define the custom theme
plot_theme <- theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 13, color = "black", margin = margin(r = 5)),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 16, face = "bold"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.margin = margin(20, 20, 20, 20),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank()
  )

# Plot top GO terms for time-based enrichment with Fold Enrichment values
ggplot(go_enrichment_time_top, aes(x = reorder(Description, FoldEnrichment), y = FoldEnrichment)) +
  geom_bar(stat = "identity", fill = "#87CEFA", color = "black", width = 0.6) +
  coord_flip() +
  labs(
    title = "Top GO Terms Enrichment (Time-based: 72H vs 24H)",
    x = "GO Term Description",
    y = "Fold Enrichment"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +  # Adds more space above bars
  plot_theme


# Plot all GO terms for condition-based enrichment
ggplot(go_enrichment_condition_df, aes(x = reorder(Description, FoldEnrichment), y = FoldEnrichment)) +
  geom_bar(stat = "identity", fill = "#FFB6C1", color = "black", width = 0.7) +
  coord_flip() +
  labs(title = "Condition-based GO Term Enrichment (SARS-CoV-2 vs Mock)",
       x = "GO Term", y = "Fold Enrichment") +
  plot_theme
