library(tidyr)
library(ggplot2)
library(viridis)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/heat_map/")

# Read the data
#data <- read.table("pacbio_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#data <- read.table("revio_gene_coverage_depth.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
data <- read.table("promethion_gene_coverage_depth.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Handle duplicate gene names
rownames(data) <- make.unique(data[, 1])  # Ensure unique gene names
gene_order <- rownames(data)              # Store the original gene order
data <- data[, -1]                        # Remove the first column (gene names)

# Transpose the data so genes are columns and samples are rows
heatmap_data <- t(as.matrix(data))        # Transpose the matrix
heatmap_data[heatmap_data > 50] <- 50
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

# Plot the heatmap with right-aligned gene labels
figure <- ggplot(df_long, aes(x = Gene, y = Sample, fill = Depth)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),  # Right-align x-axis labels
    axis.ticks.x = element_blank(),                               # Remove x-axis ticks
    axis.text.y = element_text(size = 7),                        # Adjust y-axis text size
    panel.grid = element_blank(),                                 # Remove gridlines
    legend.title = element_text(size = 10),                       # Adjust legend title size
    legend.text = element_text(size = 10)                         # Adjust legend text size
  ) +
  labs(
    x = "MHC Genes",                  # Label for x-axis
    y = "Samples",                    # Label for y-axis
    fill = "Mean\nCoverage\nDepth"    # Legend title
  ) +
  scale_x_discrete(labels = levels(df_long$Gene))  # Apply all gene labels for x-axis

figure

ggsave(filename = "new_heat_map_red.png", plot = figure, width=169, units = "mm")

#==================
library(tidyr)
library(ggplot2)
library(viridis)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/heat_map/")

revio <- read.table("revio_gene_coverage_depth.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
prom <- read.table("promethion_gene_coverage_depth.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

sample_order <- rev(colnames(revio)[-1])

# Helper to reshape one table
reshape_for_heatmap <- function(data, platform_name) {
  rownames(data) <- make.unique(data[, 1])
  gene_order <- rownames(data)
  data <- data[, -1]
  mat <- t(as.matrix(data))
  mat[mat > 50] <- 50
  df <- as.data.frame(mat)
  df$Sample <- rownames(mat)
  df_long <- pivot_longer(df, cols = -Sample, names_to = "Gene", values_to = "Depth")
  df_long$Platform <- platform_name
  df_long$Depth <- as.numeric(df_long$Depth)
  df_long$Sample <- factor(df_long$Sample, levels = sample_order)  # Preserve order here
  df_long$Gene <- factor(df_long$Gene, levels = gene_order)        # (Optional) gene order too
  df_long
}

# Combine into one long dataframe
df_long <- rbind(
  reshape_for_heatmap(revio, "Revio"),
  reshape_for_heatmap(prom, "PromethION")
)

df_long$Platform <- factor(df_long$Platform, levels = c("Revio", "PromethION"))

# Plot
figure2 <- ggplot(df_long, aes(x = Gene, y = Sample, fill = Depth)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "viridis", direction = -1) +
  theme_minimal() +
  facet_wrap(~ Platform, ncol = 1, labeller = as_labeller(c(Revio = "PacBio Revio", PromethION = "ONT PromethION"))) +
  theme(
    axis.text.x = element_text(size = 5, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 4.5),
    strip.text = element_text(size = 8),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.25, "cm"),
    legend.key.width = unit(2, "cm"),
    legend.position = "top"
  ) +
  labs(
    x = "MHC Genes",
    y = "Samples",
    fill = "Mean Coverage Depth"
  )

figure2

ggsave(filename = "combined_heat_map.png", plot = figure2)
ggsave(filename = "combined_heat_map.pdf", plot = figure2)

#==============================================================================
library(tidyr)
library(ggplot2)
library(viridis)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/heat_map/")

# Load data
revio <- read.table("revio_gene_coverage_depth.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
prom <- read.table("promethion_gene_coverage_depth.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Define desired sample order (reverse of column names, excluding first column)
sample_order <- colnames(revio)[-1]

# Helper to reshape one table
reshape_for_heatmap <- function(data, platform_name, sample_order) {
  gene_names <- make.unique(data[[1]])  # ensure unique gene names
  rownames(data) <- gene_names
  gene_order <- rev(gene_names)
  data <- data[, -1]
  mat <- t(as.matrix(data))
  mat[mat > 50] <- 50
  df <- as.data.frame(mat)
  df$Sample <- rownames(mat)
  df_long <- pivot_longer(df, cols = -Sample, names_to = "Gene", values_to = "Depth")
  df_long$Platform <- platform_name
  df_long$Depth <- as.numeric(df_long$Depth)
  df_long$Sample <- factor(df_long$Sample, levels = sample_order)  # enforce x-axis order
  df_long$Gene <- factor(df_long$Gene, levels = gene_order)        # enforce y-axis order
  df_long
}

# Reshape each dataset
revio_long <- reshape_for_heatmap(revio, "Revio", sample_order)
prom_long <- reshape_for_heatmap(prom, "PromethION", sample_order)

# Plot for Revio
plot_revio <- ggplot(revio_long, aes(x = Sample, y = Gene, fill = Depth)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#0072B2") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(size=7, face = "bold"),
    axis.title.y = element_text(size=7, face="bold"),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 7),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.position = "right"
  ) +
  labs(x = "Samples", y = "MHC Genes", fill = "Mean\nCoverage\nDepth")

# Plot for PromethION
plot_prom <- ggplot(prom_long, aes(x = Sample, y = Gene, fill = Depth)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "#0072B2") + 
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
    axis.title.x = element_text(size=7, face = "bold"),
    axis.title.y = element_text(size=7, face="bold"),
    axis.text.y = element_text(size = 5),
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 7),
    legend.key.height = unit(1.5, "cm"),
    legend.key.width = unit(0.25, "cm"),
    legend.position = "right"
  ) +
  labs(x = "Samples", y = "MHC Genes", fill = "Mean\nCoverage\nDepth")

# Display both plots
plot_prom
plot_revio

ggsave(filename = "revio_heat_map.pdf", plot = plot_revio)
ggsave(filename = "promethion_heat_map.pdf", plot = plot_prom)

#==============================================================================




