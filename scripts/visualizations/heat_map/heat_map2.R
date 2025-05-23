library(tidyr)
library(ggplot2)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/heat_map/")

# Read the data
data <- read.table("pacbio_genes_all.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Handle duplicate gene names
rownames(data) <- make.unique(data[, 1])  # Ensure unique gene names
gene_order <- rownames(data)              # Store the original gene order
data <- data[, -1]                        # Remove the first column (gene names)

# Transpose the data so genes are columns and samples are rows
heatmap_data <- t(as.matrix(data))        # Transpose the matrix
heatmap_data[heatmap_data > 100] <- 100
sample_order <- rev(rownames(heatmap_data))  # Reverse the sample order

# Convert to a data frame for ggplot
df <- as.data.frame(heatmap_data)
df$Sample <- rownames(heatmap_data)       # Add sample names as a column

# Pivot to long format for ggplot
df_long <- pivot_longer(
  df,
  cols = -Sample,                  # All columns except "Sample"
  names_to = "Gene",               # Gene names
  values_to = "Depth"              # Depth values
)

# Ensure Samples and Genes retain their original order
df_long$Sample <- factor(df_long$Sample, levels = sample_order)  # Reverse sample order
df_long$Gene <- factor(df_long$Gene, levels = gene_order)  

# Specify genes to label
genes_to_label <- c(
  "HLA-A", "HLA-B", "HLA-C", "RPL15P4", "BTNL2", "HLA-DRB1", "HLA-DPA1"
)

# Conditional axis labels for specified genes
gene_labels <- ifelse(levels(df_long$Gene) %in% genes_to_label, levels(df_long$Gene), "")

# Plot the heatmap
figure <- ggplot(df_long, aes(x = Gene, y = Sample, fill = Depth)) +
  geom_tile(color = "white") +
  #scale_fill_gradient(low = "white", high = "red") +
  scale_fill_viridis_c(option = "rocket", name = "Mean\nCoverage\nDepth") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0.5),  # Rotate x-axis labels
    axis.ticks.x = element_blank(),                               # Remove x-axis ticks
    axis.text.y = element_text(size = 10),                        # Adjust y-axis text size
    panel.grid = element_blank(),                                 # Remove gridlines
    legend.title = element_text(size = 10),                       # Adjust legend title size
    legend.text = element_text(size = 10)                         # Adjust legend text size
  ) +
  labs(
    x = "MHC Genes",                  # Label for x-axis
    y = "Samples",                    # Label for y-axis
    fill = "Mean\nCoverage\nDepth"    # Legend title
  ) +
  scale_x_discrete(labels = gene_labels)  # Apply conditional labels for x-axis

figure


ggsave(filename = "heat_map_rocket.png", plot = figure, width=169, units = "mm")
