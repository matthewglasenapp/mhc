# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork) # For stitching plots together

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/concordance/")

### First Plot ###
# Read the first CSV file
data1 <- read.csv("concordance_results.csv", header = TRUE, stringsAsFactors = FALSE)

# Clean up the numeric columns
data1$recall <- as.numeric(gsub("[^0-9.]", "", data1$recall))
data1$precision <- as.numeric(gsub("[^0-9.]", "", data1$precision))
data1$f1 <- as.numeric(gsub("[^0-9.]", "", data1$f1))

# Rename variables
data1$type[data1$type == "snp"] <- "SNV"
data1$type[data1$type == "indel"] <- "Indel"
data1$platform[data1$platform == "revio"] <- "PacBio Revio"
data1$platform[data1$platform == "promethion"] <- "ONT PromethION"

# Reorder levels
data1$type <- factor(data1$type, levels = c("SNV", "Indel"))
data1$platform <- factor(data1$platform, levels = c("PacBio Revio", "ONT PromethION"))

# Reshape the data to long format
data1_long <- data1 %>%
  pivot_longer(cols = c(recall, precision, f1),
               names_to = "metric",
               values_to = "value")

# Standardize the names of `metric`
data1_long$metric <- recode(data1_long$metric,
                            recall = "Recall",
                            precision = "Precision",
                            f1 = "F1")

# Define color and shape mappings
cb_palette <- c("Recall" = "#D55E00", "Precision" = "#0072B2", "F1" = "#009E73")
shape_mapping <- c("Recall" = 17, "Precision" = 15, "F1" = 16)

# Create the first plot
plot1 <- ggplot(data1_long, aes(x = sample, y = value, shape = metric, color = metric, group = metric)) +
  geom_point(size = 3) +
  facet_grid(type ~ platform, scales = "free_y") +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Metric Value",
    shape = "Metric",
    color = "Metric"
  ) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = cb_palette) +
  theme(
    axis.text.x = element_text(hjust = 0.5),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    strip.text = element_text(size = 8, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 8)
  )

### Second Plot ###
# Read the second CSV file
data2 <- read.csv("compare_concordance.csv", header = TRUE, stringsAsFactors = FALSE)

# Clean up the numeric columns
data2$recall <- as.numeric(gsub("[^0-9.]", "", data2$recall))
data2$precision <- as.numeric(gsub("[^0-9.]", "", data2$precision))
data2$f1 <- as.numeric(gsub("[^0-9.]", "", data2$f1))

# Rename variables
data2$platform[data2$platform == "revio"] <- "Revio"
data2$platform[data2$platform == "promethion"] <- "PromethION"

# Reshape the data to long format
data2_long <- data2 %>%
  pivot_longer(cols = c(recall, precision, f1),
               names_to = "metric",
               values_to = "value")

# Standardize the names of `metric`
data2_long$metric <- recode(data2_long$metric,
                            recall = "Recall",
                            precision = "Precision",
                            f1 = "F1")

# Ensure genes and platforms are ordered
data2_long$gene <- factor(data2_long$gene, levels = c("HLA-B", "C4A", "C4B"))
data2_long$platform <- factor(data2_long$platform, levels = c("Revio", "PromethION"))

# Create the second plot
plot2 <- ggplot(data2_long, aes(x = sample, y = value, shape = metric, color = metric, group = metric)) +
  geom_point(size = 3) +
  facet_grid(platform ~ gene, scales = "fixed") +
  theme_bw() +
  labs(
    x = "Sample",
    y = "Metric Value",
    shape = "Metric",
    color = "Metric"
  ) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = cb_palette) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(size = 8, face = "bold"),
    strip.text.y = element_text(size = 8, face = "bold"),
    strip.background = element_rect(fill = "grey90"),
    legend.title = element_text(size = 12),
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 8)
  )

### Combine Plots with Patchwork ###
combined_plot <- plot1 / plot2 + 
  plot_layout(guides = 'collect') & 
  theme(legend.position = "bottom")

# Display the combined plot
print(combined_plot)

ggsave(filename = "concordance_combined.pdf", plot = combined_plot)
ggsave(filename = "concordance_combined.png", plot = combined_plot)

plot3 <- ggplot(data1, aes(x = recall, y = precision, color = platform)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~type) +
  theme_bw() +
  scale_color_manual(
    values = c(
      "PacBio Revio" = "#5E3C99",      # magenta
      "ONT PromethION" = "#999999"     # neutral gray
    )
  ) +
  labs(
    x = "Sensitivity (Recall)",
    y = "Precision",
    color = "Platform",
  ) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    legend.position = "top",
    legend.direction = "horizontal",
  )

plot2

ggsave(filename = "concordance.pdf", plot = plot3)


