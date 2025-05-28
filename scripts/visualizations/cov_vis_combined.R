library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# Define BED + region pairs
region_configs <- list(
  list(bed = "regions.bed",        range = c(29720774, 29958570)),
  list(bed = "regions2.bed",       range = c(30000000, 30499999)),
  list(bed = "regions3.bed",       range = c(31262000, 31389999)),
  list(bed = "class2_part1.bed",   range = c(32572500, 32770000)),
  list(bed = "class2_part2.bed",   range = c(33064568, 33089696)),
  list(bed = "probe_gaps2.bed",    range = c(31519479, 32407181))  # MHC Class III
)

# Load coverage once
revio_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"
promethion_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds"
revio_all <- readRDS(revio_file)
promethion_all <- readRDS(promethion_file)

# Processing helper
process_data <- function(data, filter_range) {
  data %>%
    mutate(window = floor(base / 100) * 100 + 50) %>%
    group_by(window) %>%
    summarize(
      mean_depth = mean(mean_depth, na.rm = TRUE),
      std_depth = mean(std_depth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    rename(base = window) %>%
    filter(base >= filter_range[1], base <= filter_range[2])
}

# Main plotting function
make_plot <- function(bed_file, filter_range, use_boxes_only = FALSE) {
  annotations <- read.table(bed_file, header = FALSE, col.names = c("chr", "start", "stop", "name")) %>%
    select(start, stop, name) %>%
    mutate(midpoint = (start + stop) / 2, name = str_replace_all(name, "\\n", "\n"))
  
  revio_data <- process_data(revio_all, filter_range) %>%
    rename(revio_mean = mean_depth, revio_std = std_depth)
  promethion_data <- process_data(promethion_all, filter_range) %>%
    rename(prom_mean = mean_depth, prom_std = std_depth)
  
  combined_data <- full_join(revio_data, promethion_data, by = "base") %>%
    replace_na(list(revio_mean = 0, revio_std = 0, prom_mean = 0, prom_std = 0)) %>%
    rowwise() %>%
    mutate(
      platforms = list(if (revio_mean < prom_mean) c("PromethION", "Revio") else c("Revio", "PromethION")),
      depths = list(if (revio_mean < prom_mean) c(prom_mean, revio_mean) else c(revio_mean, prom_mean)),
      stds = list(if (revio_mean < prom_mean) c(prom_std, revio_std) else c(revio_std, prom_std))
    ) %>%
    ungroup() %>%
    select(base, platforms, depths, stds) %>%
    unnest(cols = c(platforms, depths, stds)) %>%
    rename(platform = platforms, mean_depth = depths, std_depth = stds)
  
  max_y <- max(combined_data$mean_depth)
  rect_ymin <- max_y * 1.02
  rect_ymax <- max_y * 1.07
  text_y <- max_y * 1.09
  
  coverage_plot <- ggplot(combined_data, aes(x = base, y = mean_depth, fill = platform)) +
    geom_bar(stat = "identity", position = "identity", width = 80, alpha = 1) +
    scale_fill_manual(values = c("Revio" = "blue", "PromethION" = "red")) +
    xlab("Position on Chromosome 6 (Mb)") +
    ylab("Mean\nCoverage\nDepth") +
    scale_x_continuous(
      limits = filter_range,
      labels = function(x) sprintf("%.2f", x / 1e6),
      breaks = pretty(filter_range)
    ) +
    scale_y_continuous(limits = c(0, text_y * 1.1), expand = c(0, 0)) +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.2),         # thinner tick marks
      axis.text.x = element_text(size = 4),
      axis.text.y = element_text(size = 4),
      axis.title.x = element_text(size = 4, face = "bold"),
      axis.title.y = element_text(size = 4, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),
      legend.position = "none",
      plot.margin = margin(10, 10, 6, 10)
    )
  
  if (use_boxes_only) {
    coverage_plot <- coverage_plot +
      geom_rect(data = annotations,
                aes(xmin = start, xmax = stop, ymin = 0, ymax = text_y * 1.1),
                inherit.aes = FALSE,
                fill = "grey90", color = NA)
  } else {
    coverage_plot <- coverage_plot +
      geom_rect(data = annotations,
                aes(xmin = start, xmax = stop, ymin = rect_ymin, ymax = rect_ymax),
                inherit.aes = FALSE, linewidth = 0.1, color = "black", fill = "black") +
      geom_text(data = annotations,
                aes(x = midpoint, y = text_y, label = name),
                inherit.aes = FALSE, size = 1, vjust = 0, color = "black", face = "bold")
  }
  
  return(coverage_plot)
}

# Build all plots (flag only the last one as "boxes only")
plots <- map2(
  region_configs,
  seq_along(region_configs),
  ~ make_plot(.x$bed, .x$range, use_boxes_only = .y == 6)
)

# Assemble layout: 1 + 2 + 2 + 1
final_plot <- (
  plots[[1]] /
    (plots[[2]] | plots[[3]]) /
    (plots[[4]] | plots[[5]]) /
    plots[[6]]
) + plot_layout(heights = c(1, 1, 1, 1))

# Display
print(final_plot)

# Save
ggsave(filename = "delete.pdf", plot = final_plot, width = 169, units = "mm")

