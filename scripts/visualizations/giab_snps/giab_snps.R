# Load necessary library
library(ggplot2)

# Read the data
data <- read.csv("HG002_gene_snps.csv")

# Reverse the order of the genes for the Y-axis
data$gene <- factor(data$gene, levels = rev(data$gene))

# Create the horizontal bar chart with reversed Y-axis order
figure <- ggplot(data, aes(x = snps, y = gene)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(
    x = "# SNPs (GIAB)",
    y = "Gene"
  ) +
  theme(
    axis.text.y = element_text(size = 5, margin = margin(r = -25)), # Adjust label proximity
    axis.title = element_text(size = 9, face = "bold"),
    panel.grid.major.y = element_blank(), # Remove default gridlines
    panel.grid.minor.y = element_blank()
  )

figure

# Save the plot
ggsave(filename = "giab_snps.pdf", plot = figure)
ggsave(filename = "giab_snps.png", plot = figure)


