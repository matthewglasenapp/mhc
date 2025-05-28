library(tidyverse)
library(ggplot2)
library(patchwork)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# Example BED file with chr, start, stop, and name
bed_file <- "regions.bed"
#bed_file <- "regions2.bed"
#bed_file <- "regions3.bed"
#bed_file <- "class2_part1.bed"
#bed_file <- "class2_part2.bed"
#bed_file <- "probe_gaps2.bed"
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
#filter_range <- c(32572500, 32770000)

# MHC Class II Part 2
#filter_range <- c(33064568, 33129084)

# MHC Class III
#filter_range <- c(31519479, 32407181)

# Process data for Revio and PromethION
revio_data <- process_data(revio_file, filter_range)
promethion_data <- process_data(promethion_file, filter_range)

# Add platform label and prepare for alignment
revio_data <- revio_data %>% rename(revio_mean = mean_depth, revio_std = std_depth)
promethion_data <- promethion_data %>% rename(prom_mean = mean_depth, prom_std = std_depth)

# Merge and build long format with dynamic draw order
combined_wide <- full_join(revio_data, promethion_data, by = "base") %>%
  replace_na(list(revio_mean = 0, revio_std = 0, prom_mean = 0, prom_std = 0))

combined_data <- combined_wide %>%
  rowwise() %>%
  mutate(
    platforms = list(
      if (revio_mean < prom_mean) {
        c("PromethION", "Revio")
      } else {
        c("Revio", "PromethION")
      }
    ),
    depths = list(
      if (revio_mean < prom_mean) {
        c(prom_mean, revio_mean)
      } else {
        c(revio_mean, prom_mean)
      }
    ),
    stds = list(
      if (revio_mean < prom_mean) {
        c(prom_std, revio_std)
      } else {
        c(revio_std, prom_std)
      }
    )
  ) %>%
  ungroup() %>%
  select(base, platforms, depths, stds) %>%
  unnest(cols = c(platforms, depths, stds)) %>%
  rename(platform = platforms, mean_depth = depths, std_depth = stds)

# Create combined bar plot with overlay
coverage_plot <- ggplot(combined_data, aes(x = base, y = mean_depth, fill = platform)) +
  geom_bar(stat = "identity", position = "identity", width = 80, alpha = 1) +
  scale_fill_manual(values = c("Revio" = "blue", "PromethION" = "red")) +
  xlab("Position on Chromosome 6 (Mb)") +
  ylab("Mean\nCoverage\nDepth") +
  scale_x_continuous(
    limits = filter_range,
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(filter_range, n = 10)
  ) +
  scale_y_continuous(
    limits = c(0, max(combined_data$mean_depth) * 1.1),
    expand = c(0, 0)
  ) +
  theme(
    panel.background = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.ticks = element_line(linewidth = 0.5),
    axis.title.x = element_text(size = 12, face = "bold"),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(angle = 0, size = 8, hjust = 0.5, vjust = 0.5, face = "bold"),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(0.2, "cm")
  )


# Create annotations plot
annotation_plot <- ggplot() +
  geom_rect(
    data = annotations,
    aes(xmin = start, xmax = stop, ymin = 0.5, ymax = 1),
    inherit.aes = FALSE,
    linewidth = 0.3,  # thinner border
    color = "black",
    fill = NA
  ) +
  geom_text(
    data = annotations,
    aes(x = midpoint, y = 0.25, label = name),
    inherit.aes = FALSE,
    size = 2.5,       # smaller text
    vjust = 1,
    color = "black"
  ) +
  scale_x_continuous(
    limits = filter_range,
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(filter_range, n = 10)
  ) +
  scale_y_continuous(limits = c(0, 1.5), expand = c(0, 0)) +
  theme_void() +  # remove all background axis elements
  theme(
    plot.margin = margin(2, 10, 0, 10)  # minimal top spacing
  )

# Combine annotation plot and coverage plot vertically
combined_plot <- annotation_plot / coverage_plot +
  plot_layout(heights = c(1, 5), guides = "collect") & 
  theme(
    plot.margin = margin(2, 10, 2, 10),
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 8),
    legend.spacing.x = unit(0.2, "cm")
  )


# Display the combined plot
print(combined_plot)


ggsave(filename = "delete.pdf", plot = combined_plot, width=169, units = "mm")

