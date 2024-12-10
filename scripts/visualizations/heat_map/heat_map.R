library(pheatmap)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/heat_map/")

# Read the data
data <- read.table("pacbio_genes.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Select numeric columns for the heatmap
heatmap_data <- as.matrix(data[, 2:ncol(data)])

# Transpose the matrix to put genes on the x-axis
heatmap_data <- t(heatmap_data)
heatmap_data[heatmap_data > 1000] <- 1000
rownames(heatmap_data) <- colnames(data)[2:ncol(data)]  # Set sample names as row labels

# Generate the heatmap with a legend title
pheatmap(
  heatmap_data,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "red"))(50),  # White-to-red gradient
  show_rownames = TRUE,          # Show row labels (samples) on the left
  show_colnames = FALSE,         # Suppress column labels (genes)
  na_col = "white",              # Set background for NA values explicitly
  legend_breaks = c(0, 250, 500, 1000),  # Optional: Define legend breaks
  legend_labels = c("0", "250", "500", "1000"),  # Add labels to the legend
)


# Save the heatmap to a PDF
pdf("gene_heatmap.pdf", width = 10, height = 8)  # Dimensions in inches
grid.draw(heatmap$gtable)  # Draw the heatmap using grid
dev.off()  # Close the PDF device

