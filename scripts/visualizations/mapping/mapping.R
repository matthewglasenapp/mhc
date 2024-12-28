# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(patchwork)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/mapping/")

# File paths
revio_file <- "revio_flagstat_results.csv"
promethion_file <- "promethion_flagstat_results.csv"

# Function to process and clean data
process_data <- function(file_path) {
  # Read data
  data <- read_csv(file_path, col_names = TRUE) %>%
    mutate(
      total = as.numeric(gsub(",", "", total)),  # Remove commas and convert to numeric
      duplicates = as.numeric(gsub(",", "", duplicates)),  # Remove commas and convert to numeric
      unique_reads = total - duplicates  # Calculate unique reads
    ) %>%
    select(sample, unique_reads) %>%  # Keep only necessary columns
    mutate(
      sample = factor(sample, levels = unique(sample))  # Retain the default order of samples in the data
    ) %>%
    drop_na(unique_reads)  # Remove rows with NA in unique_reads
  
  return(data)
}

# Process Revio data
revio_data <- process_data(revio_file)

# Process Promethion data and match sample order to Revio
promethion_data <- process_data(promethion_file) %>%
  mutate(sample = factor(sample, levels = levels(revio_data$sample)))  # Match sample order

# Debug: Verify processed data alignment
message("Revio Data:")
print(head(revio_data))
message("Promethion Data:")
print(head(promethion_data))

# Function to create bar plots
create_bar_plot <- function(data, color, y_label, remove_x = FALSE) {
  ggplot(data, aes(x = sample, y = unique_reads)) +
    geom_bar(stat = "identity", fill = color, color = color, alpha = 0.8) +
    scale_y_continuous(labels = scales::scientific) +  # Scientific notation for y-axis
    ylab(y_label) +
    theme_minimal() +
    theme(
      axis.text.x = if (remove_x) element_blank() else element_text(angle = 45, hjust = 1, size = 8),
      axis.title.x = if (remove_x) element_blank() else element_text(size = 12),
      axis.ticks.x = if (remove_x) element_blank() else element_line(size = 0.5),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, face = "bold"),
      panel.grid.major.x = element_blank()  # Remove vertical grid lines
    )
}

# Create plots
revio_plot <- create_bar_plot(revio_data, "blue", "Revio\nUnique\nReads", remove_x = TRUE)
promethion_plot <- create_bar_plot(promethion_data, "red", "PromethION\nUnique\nReads")

# Combine plots vertically using patchwork
combined_plot <- revio_plot / promethion_plot

# Display the combined plot
print(combined_plot)

ggsave(filename = "mapping_bar_chart.pdf", plot = combined_plot, width=169, units = "mm")
ggsave(filename = "mapping_bar_chart.png", plot = combined_plot, width=169, units = "mm")
