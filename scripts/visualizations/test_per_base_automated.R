library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# Example BED file with chr, start, stop, and name
bed_file <- "regions.bed"
#bed_file <- "probe_gaps.bed"
#bed_file <- "probe_gaps.bed"
#bed_file <- "probe_gaps.bed"
annotations <- read.table(bed_file, header = FALSE, col.names = c("chr", "start", "stop", "name"))

# Process annotation data
annotations <- annotations %>%
  select(start, stop, name) %>%
  mutate(
    midpoint = (start + stop) / 2,
    name = str_replace_all(name, "\\\\n", "\n")  # Replace \n with actual newlines
  )

# File paths for per-base coverage files
revio_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"
promethion_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds"

# Function to process and filter data
process_data <- function(file_path, filter_range) {
  data <- readRDS(file_path)
  collapsed_data <- data %>%
    mutate(window = floor(base / 100) * 100 + 50) %>%
    group_by(window) %>%
    summarize(
      mean_depth = mean(mean_depth, na.rm = TRUE),
      std_depth = mean(std_depth, na.rm = TRUE)
    ) %>%
    rename(base = window) %>%
    filter(base >= filter_range[1] & base <= filter_range[2])
  
  return(collapsed_data)
}

# Filter range for visualization
# Modify this to toggle between different regions of interest
# MHC Class I Part 1
filter_range <- c(29720774, 30016539)

# MHC Class I Part 2
#filter_range <- c(30000000, 30499999)

# MHC Class I Part 3
#filter_range <- c(31262000, 31389999)

# MHC Class II Part 1
# filter_range <- c(32572500, 32770000)

# MHC Class II Part 2
# filter_range <- c(33064568, 33129084)

# MHC Class III
#filter_range <- c(31519479, 32407181)


# Process data for Revio and Promethion
revio_data <- process_data(revio_file, filter_range)
promethion_data <- process_data(promethion_file, filter_range)

# Function to create bar plots
create_bar_plot <- function(data, color, y_label, remove_x = FALSE, remove_x_line = FALSE) {
  plot <- ggplot(data, aes(x = base, y = mean_depth)) +
    geom_bar(stat = "identity", fill = color, color = color, width = 80, alpha = 0.8) +
    xlab("Position on Chromosome 6 (Mb)") +
    ylab(y_label) +
    theme(
      panel.background = element_blank(),
      panel.border = element_blank(),
      axis.line.y = element_line(colour = "black", size = 1),
      axis.line.x = if (remove_x_line) element_blank() else element_line(colour = "black", size = 1),
      axis.text.x = if (remove_x) element_blank() else element_text(size = 12, angle = 45, hjust = 1),
      axis.ticks.x = if (remove_x) element_blank() else element_line(size = 0.5),
      axis.title.x = if (remove_x) element_blank() else element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 8),
      axis.title.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, face = "bold")
    ) +
    scale_x_continuous(
      limits = filter_range,
      labels = function(x) sprintf("%.2f", x / 1e6),
      breaks = pretty(filter_range, n = 10)
    ) +
    scale_y_continuous(
      limits = c(0, max(data$mean_depth) * 1.1),
      expand = c(0, 0)
    )
  
  return(plot)
}

# Create bar plots for Revio and Promethion
revio_plot <- create_bar_plot(revio_data, "blue", "Revio\nMean\nCoverage\nDepth", remove_x = TRUE, remove_x_line = TRUE)
promethion_plot <- create_bar_plot(promethion_data, "red", "PromethION\nMean\nCoverage\nDepth")

# Create annotations plot
annotation_plot <- ggplot() +
  geom_rect(
    data = annotations,
    aes(xmin = start, xmax = stop, ymin = 10, ymax = 13),
    inherit.aes = FALSE,
    size = 0.5, color = "black", fill = NA
  ) +
  geom_text(
    data = annotations,
    aes(x = midpoint, y = 5, label = name),
    inherit.aes = FALSE,
    color = "black", vjust = 0, size = 3
  ) +
  scale_x_continuous(
    limits = filter_range,
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(filter_range, n = 10)
  ) +
  scale_y_continuous(
    limits = c(0, 20),
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
combined_plot <- revio_plot / annotation_plot / promethion_plot

# Display the combined plot
print(combined_plot)


ggsave(filename = "HLA_III.pdf", plot = combined_plot, width=169, units = "mm")
