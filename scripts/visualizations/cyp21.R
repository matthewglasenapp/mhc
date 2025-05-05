library(tidyverse)
library(ggplot2)
library(gggenes)
library(patchwork)

# Set working directory
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# File paths
revio_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"
promethion_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds"
bed_file <- "class_III_rccs.bed"

# New filter range: 32,004,688 – 32,042,644
filter_range <- c(32004688, 32042644)

# Target genes
target_genes <- c("CYP21A1P", "TNXA", "STK19B", "C4B", "CYP21A2")

# Read and filter gene annotations
annotations <- read.table(
  bed_file, header = TRUE, sep = "\t",
  col.names = c("chr", "start", "stop", "gene", "strand", "y")
) %>%
  filter(gene %in% target_genes) %>%
  filter(start <= filter_range[2] & stop >= filter_range[1]) %>%
  arrange(start) %>%
  mutate(
    gene = factor(gene, levels = target_genes),
    y = 1
  )

# Function to process coverage data
process_data <- function(file_path, filter_range, platform_name) {
  readRDS(file_path) %>%
    mutate(window = floor(base / 100) * 100 + 50) %>%
    group_by(window) %>%
    summarize(
      mean_depth = mean(mean_depth, na.rm = TRUE),
      std_depth = mean(std_depth, na.rm = TRUE)
    ) %>%
    rename(base = window) %>%
    filter(base >= filter_range[1] & base <= filter_range[2]) %>%
    mutate(platform = platform_name)
}

# Process coverage
revio_data <- process_data(revio_file, filter_range, "Revio")
promethion_data <- process_data(promethion_file, filter_range, "PromethION")

# Combine and format for plotting
combined_data <- bind_rows(revio_data, promethion_data) %>%
  mutate(platform = factor(platform, levels = c("PromethION", "Revio")))

# Coverage plot
coverage_plot <- ggplot(combined_data, aes(x = base, y = mean_depth, fill = platform)) +
  geom_bar(stat = "identity", position = "stack", width = 100, alpha = 1) +
  scale_fill_manual(values = c("PromethION" = "blue", "Revio" = "red")) +
  xlab("Position on Chromosome 6 (Mb)") +
  ylab("Mean Coverage Depth") +
  scale_x_continuous(
    limits = filter_range,
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(filter_range, n = 6)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.title = element_blank(),
    legend.text = element_text(size = 9)
  )

# Gene annotation plot (gggenes)
annotation_plot <- ggplot(annotations, aes(
  xmin = start, xmax = stop, y = y, fill = gene, forward = strand == "+"
)) +
  geom_hline(yintercept = 1, color = "grey", linetype = "solid", size = 0.5) +
  geom_gene_arrow(color = "black", size = 0.7, alpha = 1) +
  geom_gene_label(aes(label = gene), size = 4) +
  scale_fill_brewer(palette = "Set3") +
  scale_x_continuous(limits = filter_range) +  # Match coverage x-axis
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.6, "cm"),
    legend.margin = margin(10, 10, 10, 10),
    plot.margin = margin(10, 10, 10, 10)
  )

# Combine plots
final_plot <- annotation_plot / coverage_plot +
  plot_layout(heights = c(1, 2), guides = "collect")

# Display
print(final_plot)


ggsave(filename = "cyp.pdf", plot = final_plot)

