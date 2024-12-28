library(tidyverse)
library(ggplot2)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/")

# Example BED file with chr, start, stop, and name
bed_file <- "regions3.bed"
annotations <- read.table(bed_file, header = FALSE, col.names = c("chr", "start", "stop", "name"))

# Ignore the first column (chr) and work with start, stop, and name
annotations <- annotations %>%
  select(start, stop, name) %>%
  mutate(midpoint = (start + stop) / 2)

# per-base coverage file
hla_per_base <- "/Users/matt/Documents/GitHub/mhc/clean_data/revio_mean_std_depth.rds"

data <- readRDS(hla_per_base)

# Collapse data into 100-base-pair windows
collapsed_data <- data %>%
  mutate(window = floor(base / 100) * 100 + 50) %>% # Create windows centered at midpoints
  group_by(window) %>%
  summarize(
    mean_depth = mean(mean_depth),
    std_depth = mean(std_depth)
  )

# Rename 'window' to 'base' for consistency
collapsed_data <- collapsed_data %>%
  rename(base = window)

# MHC Class I Part 1
#data_filtered <- collapsed_data %>% filter(base >= 29720774 & base <= 30016539)
# MHC Class I Part 2
#data_filtered <- collapsed_data %>% filter(base >= 30000000 & base <= 30499999)
# MHC Class I Part 3
data_filtered <- collapsed_data %>% filter(base >= 31262000 & base <= 31389999)
#MHC Class II Part 1
#data_filtered <- collapsed_data %>% filter(base >= 32572500 & base <= 32770000)
# MHC Class II Part 2
#data_filtered <- collapsed_data %>% filter(base >= 33064568 & base <= 33129084)
# MHC Class III
#data_filtered <- collapsed_data %>% filter(base >= 31519479 & base <= 32407181)
 
# Create the plot
figure <- ggplot(data_filtered, aes(x = base, y = mean_depth)) +
  geom_line(size = 0.75) +
  xlab("Position on Chromosome 6 (Mb)") +
  ylab("Mean Coverage Depth") +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank(),
        axis.line = element_line(colour = "black", linewidth = 1),
        axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 12, face = "bold")) +
  scale_x_continuous(
    labels = function(x) sprintf("%.2f", x / 1e6),
    breaks = pretty(data_filtered$base, n = 10)
  ) +
  scale_y_continuous(
    limits = c(-110, max(data_filtered$mean_depth)),  # Adjust limits if necessary
    expand = c(0, 0)
  ) + 
  geom_rect(
    data = annotations,
    aes(xmin = start, xmax = stop, ymin = -50, ymax = -25),
    inherit.aes = FALSE, # Prevent inheritance of aesthetics from ggplot()
    size = 0.5, color = "black", fill = NA
  ) +
  geom_text(
    data = annotations,
    aes(x = midpoint, y = -100, label = name),
    inherit.aes = FALSE, # Prevent inheritance of aesthetics from ggplot()
    color = "black", vjust = 0, size = 3
  )

# Display the figure
print(figure)

ggsave(filename = "class_II_part1.png", plot = figure, width=169, units = "mm")

