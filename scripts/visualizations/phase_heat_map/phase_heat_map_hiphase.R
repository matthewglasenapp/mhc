library(tidyr)
library(ggplot2)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/phase_heat_map/")

# Load the data
data <- read.csv("phase_map_hiphase.csv", header = TRUE, stringsAsFactors = FALSE)

# Replace dots with dashes in column names
colnames(data) <- gsub("\\.", "-", colnames(data))

df_long <- pivot_longer(
  data,
  cols = -sample,       # Columns to pivot: All except "sample"
  names_to = "Gene",    # New column for gene names
  values_to = "Value"   # New column for gene values (1 or 0)
)

# Ensure samples and genes retain their order
df_long$sample <- factor(df_long$sample, levels = rev(unique(data$sample)))  # Reverse sample order
df_long$Gene <- factor(df_long$Gene, levels = colnames(data)[-1])  # Gene order

# Plot the heatmap
figure <- ggplot(df_long, aes(x = Gene, y = sample, fill = Value)) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "red", na.value = "white") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "none"  # Remove the legend
  ) +
  labs(
    x = "Gene",
    y = "Sample"
  )

# Print the figure
print(figure)

# Save the figure
ggsave(filename = "hiphase_heat_map.pdf", plot = figure)
ggsave(filename = "hiphase_heat_map.png", plot = figure)



