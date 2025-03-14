# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/concordance/")

# Read the CSV file
data <- read.csv("compare_concordance.csv", header = TRUE, stringsAsFactors = FALSE)

# Clean up the numeric columns (if necessary)
data$recall <- as.numeric(gsub("[^0-9.]", "", data$recall))
data$precision <- as.numeric(gsub("[^0-9.]", "", data$precision))
data$f1 <- as.numeric(gsub("[^0-9.]", "", data$f1))

# Rename 'platform' variables directly
data$platform[data$platform == "revio"] <- "Revio"
data$platform[data$platform == "promethion"] <- "PromethION"

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = c(recall, precision, f1),
               names_to = "metric",
               values_to = "value")

# Standardize the names of `metric`
data_long$metric <- recode(data_long$metric,
                           recall = "Recall",
                           precision = "Precision",
                           f1 = "F1")

# Ensure genes are ordered as columns (left to right: HLA-B, C4A, C4B)
data_long$gene <- factor(data_long$gene, levels = c("HLA-B", "C4A", "C4B"))

# Ensure platforms are ordered as rows (top to bottom: Revio, PromethION)
data_long$platform <- factor(data_long$platform, levels = c("Revio", "PromethION"))

# Debugging: Check the number of points per sample per gene per platform
debug_summary <- data_long %>%
  group_by(gene, sample, platform) %>%
  summarise(num_points = n(), .groups = "drop")

# Output debugging information
print("Debugging Summary (should be 3 data points per sample per platform):")
print(debug_summary)

# Write the debugging summary to a CSV file for further inspection if needed
write.csv(debug_summary, "debug_summary.csv", row.names = FALSE)

# Define a color-blind-friendly palette
cb_palette <- c("Recall" = "#D55E00", "Precision" = "#0072B2", "F1" = "#009E73")

# Map specific shapes to each metric
shape_mapping <- c("Recall" = 17, "Precision" = 15, "F1" = 16) # Triangle, square, circle

# Plot the data
figure <- ggplot(data_long, aes(x = sample, y = value, shape = metric, color = metric, group = metric)) +
  geom_point(size = 3) +
  facet_grid(platform ~ gene, scales = "fixed") + # Arrange genes as columns, platforms as rows, fixed scales
  theme_bw() +
  labs(
    x = "Sample",
    y = "Metric Value",
    shape = "Metric",
    color = "Metric"
  ) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = cb_palette) + # Apply color-blind-friendly palette
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 10, face = "bold"),       # Bold text for gene facet labels
    strip.text.y = element_text(size = 10, face = "bold"),       # Bold text for platform facet labels
    strip.background = element_rect(fill = "grey90"), # Light background for facet labels
    legend.title = element_text(size = 12),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )

# Display the plot
print(figure)

ggsave(filename = "concordance_comparison.pdf", plot = figure, width = 169, units = "mm")
ggsave(filename = "concordance_comparison.png", plot = figure, width = 169, units = "mm")
