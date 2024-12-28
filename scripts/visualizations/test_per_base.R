library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# File paths
revio_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"
promethion_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds"

# Function to process and filter data
process_data <- function(file_path, filter_range) {
  data <- readRDS(file_path)
  collapsed_data <- data %>%
    mutate(window = floor(base / 100) * 100 + 50) %>%  # Create 100-base-pair windows centered at midpoints
    group_by(window) %>%
    summarize(mean_depth = mean(mean_depth), std_depth = mean(std_depth)) %>%
    rename(base = window)
  
  # Filter by specified range
  collapsed_data %>% filter(base >= filter_range[1] & base <= filter_range[2])
}

# Filter ranges for visualization
filter_range <- c(29715000, 30019999)

# Process data
revio_data <- process_data(revio_file, filter_range)
promethion_data <- process_data(promethion_file, filter_range)

# Function to create bar plots
create_bar_plot <- function(data, color, y_label, remove_x = FALSE, remove_x_line = FALSE) {
  plot <- ggplot(data, aes(x = base, y = mean_depth)) +
    geom_bar(stat = "identity", fill = color, color = color, width = 80, alpha = 0.8) +  # Bar chart
    xlab("Position on Chromosome 6 (Mb)") +
    ylab(y_label) +
    theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black", size = 1), # Ensure y-axis starts at 0
      axis.line.x = if (remove_x_line) element_blank() else element_line(colour = "black", size = 1),
      axis.text.x = if (remove_x) element_blank() else element_text(size = 12, angle = 45, hjust = 1),
      axis.ticks.x = if (remove_x) element_blank() else element_line(size = 0.5),
      axis.title.x = if (remove_x) element_blank() else element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, face = "bold")
    ) +
    scale_x_continuous(
      limits = filter_range, # Ensures the x-axis matches the filter range
      labels = function(x) sprintf("%.2f", x / 1e6),
      breaks = pretty(filter_range, n = 10)
    ) +
    scale_y_continuous(
      limits = if (color == "blue") c(0, max(data$mean_depth) * 1.1) else c(-220, max(data$mean_depth) * 1.1),
      expand = c(0, 0)
    )
  
  return(plot)
}

# Create bar plots
revio_plot <- create_bar_plot(revio_data, "blue", "Revio\nMean\nCoverage\nDepth", remove_x = TRUE, remove_x_line = TRUE)
promethion_plot <- create_bar_plot(promethion_data, "red", "PromethION\nMean\nCoverage\nDepth")

# Create a new plot for rectangle annotations
annotation_plot <- ggplot() +
  geom_rect(
    aes(xmin = 29792233, ymin = 10, xmax = 29793136, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29792233, 29793136)), y = 5,
    label = "HLA-V", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29941259, ymin = 10, xmax = 29949572, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29941259, 29949572)), y = 5,
    label = "HLA-A", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29722774, ymin = 10, xmax = 29738528, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29722774, 29738528)), y = 5,
    label = "HLA-F", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29800414, ymin = 10, xmax = 29802425, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29800414, 29802425)), y = 5,
    label = "HLA-P", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29826966, ymin = 10, xmax = 29831125, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29826966, 29831125)), y = 5,
    label = "HLA-G", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29887751, ymin = 10, xmax = 29890482, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29887751, 29890482)), y = 5,
    label = "HLA-H", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29926458, ymin = 10, xmax = 29929232, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29926458, 29929232)), y = 5,
    label = "HLA-K", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 29956595, ymin = 10, xmax = 29958570, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(29956595, 29958570)), y = 5,
    label = "HLA-W", color = "black", vjust = 0, size = 3
  ) +
  geom_rect(
    aes(xmin = 30006605, ymin = 10, xmax = 30009539, ymax = 13),
    size = 0.5, color = "black", fill = NA
  ) +
  annotate(
    "text", x = mean(c(30006605, 30009539)), y = 5,
    label = "HLA-J", color = "black", vjust = 0, size = 3
  ) +
  scale_x_continuous(
    limits = filter_range,
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(filter_range, n = 10)
  ) +
  scale_y_continuous(
    limits = c(0, 20),  # Set y-axis limits to focus on the rectangles
    expand = c(0, 0)
  ) +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

# Combine the plots vertically using patchwork
combined_plot <- annotation_plot / revio_plot / promethion_plot

# Display the combined plot
print(combined_plot)


ggsave(filename = "HLA_Class_I.pdf", plot = combined_plot, width=169, units = "mm")
