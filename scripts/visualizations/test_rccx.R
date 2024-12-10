library(tidyverse)
library(ggplot2)

# Set working directory (adjust as needed)
setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# Example BED file with chr, start, stop, and name
bed_file <- "class_III_rccs.bed"
annotations <- read.table(bed_file, header = FALSE, col.names = c("chr", "start", "stop", "name"))

# Ignore the first column (chr) and work with start, stop, and name
annotations <- annotations %>%
  select(start, stop, name) %>%
  mutate(
    midpoint = (start + stop) / 2,  # Calculate midpoint for x-axis
    name = gsub("\\\\n", "\n", name)  # Convert literal \n to actual newline
  )

# Assign y-levels to stagger rectangles and calculate midpoint_y
annotations <- annotations %>%
  mutate(
    ymin = seq(-40, by = -20, length.out = n()),  # Start ymin at -40 and decrement by 20
    ymax = ymin + 15,  # Set ymax based on ymin
    midpoint_y = (ymin + ymax) / 2  # Calculate the midpoint for y-axis
  )

# Per-base coverage file
hla_per_base <- "/Users/matt/Documents/GitHub/mhc/clean_data/pacbio_mean_std_depth.rds"
data <- readRDS(hla_per_base)

# Collapse data into 100-base-pair windows
collapsed_data <- data %>%
  mutate(window = floor(base / 100) * 100 + 50) %>%  # Create windows centered at midpoints
  group_by(window) %>%
  summarize(
    mean_depth = mean(mean_depth),
    std_depth = mean(std_depth)
  )

# Rename 'window' to 'base' for consistency
collapsed_data <- collapsed_data %>%
  rename(base = window)

# Filter data for the region of interest
data_filtered <- collapsed_data %>%
  filter(base >= 31970000 & base <= 32118000)

# Create the plot
figure <- ggplot(data_filtered, aes(x = base, y = mean_depth)) +
  geom_line(size = 0.75) +
  xlab("Position on Chromosome 6 (Mb)") +
  ylab("Mean Coverage Depth") +
  theme(panel.background = element_blank()) +
  theme(
    panel.border = element_blank(),
    axis.line = element_line(colour = "black", linewidth = 1),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12, face = "bold")
  ) +
  scale_x_continuous(
    labels = function(x) sprintf("%.2f", x / 1e6),  # Convert to Mb
    breaks = pretty(data_filtered$base, n = 10)
  ) +
  scale_y_continuous(
    limits = c(-200, max(data_filtered$mean_depth)),  # Extend the y-axis below 0
    breaks = c(0, pretty(seq(0, max(data_filtered$mean_depth), length.out = 5))),  # Custom y-axis breaks
    expand = c(0, 0)
  ) +
  geom_rect(
    data = annotations,
    aes(xmin = start, xmax = stop, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,  # Prevent inheritance of aesthetics from ggplot()
    size = 0.5, color = "black", fill = NA
  ) +
  geom_text(
    data = annotations,
    aes(x = stop + 1000, y = midpoint_y, label = name),  # Use midpoint_y for vertical centering
    inherit.aes = FALSE,  # Prevent inheritance of aesthetics from ggplot()
    color = "black", hjust = 0, vjust = 0.5, size = 3.25
  )

# Display the figure
print(figure)


ggsave(filename = "rccx.png", plot = figure, width=169, units = "mm")
