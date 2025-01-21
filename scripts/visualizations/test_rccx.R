library(tidyverse)
library(ggplot2)
library(gggenes)
library(patchwork)

# Set working directory (adjust as needed)
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# File paths for coverage data
revio_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"
promethion_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds"

# Filter range for visualization
filter_range <- c(31970000, 32118000)  # Region of interest

# Load gene annotations
bed_file <- "class_III_rccs.bed"
annotations <- read.table(
  bed_file, header = TRUE, sep = "\t",
  col.names = c("chr", "start", "stop", "gene", "strand", "y")
)

# Debug: Check the annotations data
print("Debugging Annotations Data:")
print(head(annotations))

# Reorder genes for plotting and preserve legend order
annotations <- annotations %>%
  mutate(
    plot_order = ifelse(gene == "STK19", 1, 0),  # STK19 plotted last
    gene = factor(gene, levels = unique(gene))  # Preserve legend order (left-to-right appearance)
  ) %>%
  arrange(plot_order, start) %>%
  mutate(y = 1)  # Keep all genes on the same y-level

# Function to process and filter coverage data
process_data <- function(file_path, filter_range, platform_name) {
  data <- readRDS(file_path) %>%
    mutate(window = floor(base / 100) * 100 + 50) %>%
    group_by(window) %>%
    summarize(
      mean_depth = mean(mean_depth, na.rm = TRUE),
      std_depth = mean(std_depth, na.rm = TRUE)
    ) %>%
    rename(base = window) %>%
    filter(base >= filter_range[1] & base <= filter_range[2]) %>%
    mutate(platform = platform_name)  # Add platform identifier
  
  return(data)
}

# Process coverage data for Revio and PromethION
revio_data <- process_data(revio_file, filter_range, "Revio")
promethion_data <- process_data(promethion_file, filter_range, "PromethION")

# Combine coverage data for stacked bar plot
combined_data <- bind_rows(revio_data, promethion_data)

# Debug: Check combined coverage data
print("Debugging Combined Coverage Data:")
print(head(combined_data))

# Explicitly reorder platforms for legend consistency
combined_data$platform <- factor(combined_data$platform, levels = c("PromethION", "Revio"))

# Create the stacked bar plot for mean coverage depth
coverage_plot <- ggplot(combined_data, aes(x = base, y = mean_depth, fill = platform)) +
  geom_bar(stat = "identity", position = "stack", width = 100, alpha = 1) +  # Align width to 100bp bins, remove transparency
  scale_fill_manual(
    values = c("PromethION" = "blue", "Revio" = "red"),  # Corrected color mapping
    labels = c("PromethION", "Revio")
  ) +
  xlab("Position on Chromosome 6 (Mb)") +
  ylab("Mean Coverage Depth") +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),  # Ensure no grid lines
    legend.position = "bottom",  # Place the legend below the plot
    legend.direction = "vertical",  # Stack legend items vertically
    legend.title = element_blank(),
    legend.text = element_text(size = 9)  # Slightly smaller font for legend
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.2f", x / 1e6),  # Convert to Mb
    breaks = pretty(filter_range, n = 10)
  ) +
  scale_y_continuous(expand = c(0, 0))

# Debug: Check the coverage plot
print("Rendering Coverage Depth Plot:")
print(coverage_plot)

# Create the gggenes annotation plot
annotation_plot <- ggplot(annotations, aes(
  xmin = start, xmax = stop, y = y, fill = gene, forward = strand == "+"
)) +
  geom_hline(yintercept = 1, color = "grey", linetype = "solid", size = 0.5) +  # Add the grey line first
  geom_gene_arrow(
    color = "black",  # Darker outline
    size = 0.7,  # Thicker line weight
    alpha = 1
  ) +
  geom_gene_label(aes(label = gene), size = 4) +  # Adjusted size for visibility
  scale_fill_brewer(palette = "Set3") +  # Use Set3 color palette
  theme_void() +  # Completely remove all axis lines and ticks
  theme(
    legend.position = "right",  # Legend on the right
    legend.title = element_blank(),
    legend.text = element_text(size = 9),  # Smaller font for gene legend
    legend.key.size = unit(0.6, "cm"),  # Slightly smaller legend boxes
    legend.margin = margin(10, 10, 10, 10),  # Add margin around legend
    plot.margin = margin(10, 10, 10, 10)  # Add margin around the plot
  )

# Debug: Check the annotation plot
print("Rendering gggenes Annotation Plot:")
print(annotation_plot)

# Combine the annotation and coverage plots
final_plot <- annotation_plot / coverage_plot +
  plot_layout(heights = c(1, 2), guides = "collect")  # Ensure guides are properly arranged

# Display the final combined plot
print(final_plot)




ggsave(filename = "rccx.png", plot = final_plot, width=169, units = "mm")
ggsave(filename = "rccx.pdf", plot = final_plot, width=169, units = "mm")
