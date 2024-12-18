# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/concordance/")

# Read the CSV file
data <- read.csv("concordance_results.csv", header = TRUE, stringsAsFactors = FALSE)
head(data)

# Clean up the numeric columns (if necessary)
data$recall <- as.numeric(gsub("[^0-9.]", "", data$recall))
data$precision <- as.numeric(gsub("[^0-9.]", "", data$precision))
data$f1 <- as.numeric(gsub("[^0-9.]", "", data$f1))

# Rename 'type' and 'platform' variables directly
data$type[data$type == "snp"] <- "SNV"
data$type[data$type == "indel"] <- "Indel"

data$platform[data$platform == "revio"] <- "PacBio Revio"
data$platform[data$platform == "promethion"] <- "ONT PromethION"

# Reorder the levels of `type` and `platform`
data$type <- factor(data$type, levels = c("SNV", "Indel")) # SNV on top
data$platform <- factor(data$platform, levels = c("PacBio Revio", "ONT PromethION")) # PacBio Revio on left

# Reshape the data to long format
data_long <- data %>%
  pivot_longer(cols = c(recall, precision, f1),
               names_to = "metric",
               values_to = "value")

# Standardize the names of `metric`
data_long$metric <- recode(data_long$metric,
                           recall = "Recall",
                           precision = "Precision",
                           f1 = "F1")  # Ensure consistent naming

# Define a color-blind-friendly palette
cb_palette <- c("Recall" = "#D55E00", "Precision" = "#0072B2", "F1" = "#009E73")

# Map specific shapes to each metric
shape_mapping <- c("Recall" = 17, "Precision" = 15, "F1" = 16) # Triangle, square, and circle

# Plot the data
figure <- ggplot(data_long, aes(x = sample, y = value, shape = metric, color = metric, group = metric)) +
  geom_point(size = 3) +
  facet_grid(type ~ platform, scales = "free_y") + # Facet using reordered variables
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
    axis.text.x = element_text(hjust = 0.5),
    strip.text = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )

# Display the plot
figure

ggsave(filename = "concordance.pdf", plot = figure, width = 169, units = "mm")
ggsave(filename = "concordance.png", plot = figure, width = 169, units = "mm")

