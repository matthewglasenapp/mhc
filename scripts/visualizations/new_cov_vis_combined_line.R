library(tidyverse)
library(ggplot2)
library(patchwork)
library(glue)  # for debug messages

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# Define BED + region pairs
region_configs <- list(
  list(bed = "regions.bed",        range = c(29719000, 29961000)),
  list(bed = "regions2.bed",       range = c(30000000, 30499999)),
  list(bed = "regions3.bed",       range = c(31262000, 31389999)),
  list(bed = "class2_part1.bed",   range = c(32572500, 32770000)),
  list(bed = "class2_part2.bed",   range = c(33064568, 33090000)),
  list(bed = "probe_gaps.bed",     range = c(31519479, 32420000))  # MHC Class III
)

# Plot titles
region_titles <- c(
  "MHC Class I", "MHC Class I (cont.)", "MHC Class I (cont.)",
  "MHC Class II", "MHC Class II (cont.)",
  "MHC Class III"
)

# Load coverage
revio_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"
promethion_file <- "/Users/matt/Documents/GitHub/mhc/clean_data/promethion_mean_std_depth.rds"
revio_all <- readRDS(revio_file)
promethion_all <- readRDS(promethion_file)

# Helper to bin and summarize
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

# Debugging function
debug_annotation_limits <- function(idx, annotations, y_limit) {
  out_of_bounds <- annotations %>%
    filter(rect_ymin > y_limit | rect_ymax > y_limit | text_y > y_limit)
  if (nrow(out_of_bounds) > 0) {
    message(glue("Plot {idx}: {nrow(out_of_bounds)} annotation(s) outside y_limit {y_limit}"))
    print(out_of_bounds)
  } else {
    message(glue("Plot {idx}: all annotations within y_limit {y_limit}"))
  }
}

# Plotting function
make_plot <- function(bed_file, filter_range, use_boxes_only = FALSE, title_text = NULL, is_plot5 = FALSE, idx = NA) {
  annotations <- read.table(bed_file, header = FALSE, col.names = c("chr", "start", "stop", "name")) %>%
    select(start, stop, name) %>%
    mutate(midpoint = (start + stop) / 2, name = str_replace_all(name, "\\n", "\n"))
  
  revio_data <- process_data(revio_all, filter_range) %>%
    rename(revio_mean = mean_depth, revio_std = std_depth)
  promethion_data <- process_data(promethion_all, filter_range) %>%
    rename(prom_mean = mean_depth, prom_std = std_depth)
  
  combined_data <- full_join(revio_data, promethion_data, by = "base") %>%
    replace_na(list(revio_mean = 0, revio_std = 0, prom_mean = 0, prom_std = 0)) %>%
    pivot_longer(
      cols = c(revio_mean, prom_mean),
      names_to = "platform",
      values_to = "mean_depth"
    ) %>%
    mutate(platform = ifelse(platform == "revio_mean", "Revio", "PromethION"))
  
  raw_max_y <- max(combined_data$mean_depth)
  capped <- use_boxes_only
  
  if (capped) {
    combined_data <- combined_data %>% mutate(mean_depth = pmin(mean_depth, 350))
  }
  
  plot_max_y <- if (capped) 350 else raw_max_y
  rect_height <- plot_max_y * 0.05
  rect_ymin_base <- plot_max_y * 1.02
  rect_ymax_base <- rect_ymin_base + rect_height
  text_y_base <- rect_ymax_base + (plot_max_y * 0.01)
  
  if (!use_boxes_only) {
    if (is_plot5) {
      stagger_height <- rect_height * 1.5
      annotations <- annotations %>%
        mutate(
          offset = (row_number() - 1) * stagger_height + 2,
          rect_ymin = rect_ymin_base + offset,
          rect_ymax = rect_ymin + rect_height,
          text_y = text_y_base + offset
        )
    } else {
      annotations <- annotations %>%
        mutate(
          rect_ymin = rect_ymin_base,
          rect_ymax = rect_ymax_base,
          text_y = text_y_base
        )
    }
    # Updated y_limit now that annotations are finalized
    y_limit <- max(max(combined_data$mean_depth), annotations$rect_ymax, annotations$text_y) * 1.05
    debug_annotation_limits(idx, annotations, y_limit)
  } else {
    y_limit <- text_y_base + (plot_max_y * 0.05)
  }
  
  coverage_plot <- ggplot(combined_data, aes(x = base, y = mean_depth, color = platform, fill = platform))
  
  # Grey boxes behind everything
  if (use_boxes_only) {
    coverage_plot <- coverage_plot +
      geom_rect(data = annotations,
                aes(xmin = start, xmax = stop, ymin = 0, ymax = text_y_base),
                inherit.aes = FALSE,
                fill = "grey90", color = NA)
  }
  
  # Main coverage layers
  coverage_plot <- coverage_plot +
    geom_area(alpha = 0.2, position = "identity") +
    geom_line(linewidth = 0.05) +
    scale_color_manual(values = c("Revio" = "#0072B2", "PromethION" = "#D55E00")) +
    scale_fill_manual(values = c("Revio" = "#0072B2", "PromethION" = "#D55E00")) +
    xlab("Position on Chromosome 6 (Mb)") +
    ylab("Mean\nCoverage\nDepth") +
    scale_x_continuous(
      limits = filter_range,
      labels = function(x) sprintf("%.2f", x / 1e6),
      breaks = pretty(filter_range)
    ) +
    scale_y_continuous(limits = c(0, y_limit), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    theme(
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.3),
      axis.ticks = element_line(linewidth = 0.2),
      axis.text.x = element_text(size = 4),
      axis.text.y = element_text(size = 4),
      axis.title.x = element_text(size = 4, face = "bold"),
      axis.title.y = element_text(size = 4, face = "bold", angle = 0, hjust = 0.5, vjust = 0.5),
      legend.position = "none",
      plot.margin = margin(2, 6, 2, 6)
    )
  
  if (!use_boxes_only) {
    coverage_plot <- coverage_plot +
      geom_rect(data = annotations,
                aes(xmin = start, xmax = stop, ymin = rect_ymin, ymax = rect_ymax),
                inherit.aes = FALSE, linewidth = 0.1, color = "black", fill = "black") +
      geom_text(data = annotations,
                aes(x = midpoint, y = text_y, label = name),
                inherit.aes = FALSE,
                size = 1.5, vjust = 0, color = "black")
  }
  
  if (!is.null(title_text)) {
    coverage_plot <- coverage_plot +
      ggtitle(title_text) +
      theme(plot.title = element_text(size = 5, hjust = 0.5))
  }
  
  return(coverage_plot)
}

# Build plots with titles
plots <- pmap(
  list(region_configs, seq_along(region_configs), region_titles),
  function(config, idx, title) {
    make_plot(config$bed, config$range,
              use_boxes_only = idx == 6,
              title_text = title,
              is_plot5 = idx == 5,
              idx = idx)
  }
)

# Layout: 1 + 2 + 2 + 1
final_plot <- (
  plots[[1]] /
    (plots[[2]] | plots[[3]]) /
    (plots[[4]] | plots[[5]]) /
    plots[[6]]
) + plot_layout(heights = c(1, 1, 1, 1))

# Display
print(final_plot)

# Save
ggsave(filename = "coverage_plot.pdf", plot = final_plot, dpi = 600)

ggsave("coverage_plot.png", plot = final_plot, dpi = 600, width = 6.5, units = "in")


