library(tidyverse)
library(ggplot2)

setwd("/Users/matt/Documents/GitHub/mhc/scripts/visualizations/c4/")


bed_data <- read_tsv("subset_exons.bed", col_names = c("chrom", "start", "end"))

# per-base coverage file
hla_per_base <- "/Users/matt/Documents/GitHub/mhc/clean_data/pacbio_mean_std_depth.rds"

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

data_filtered <- collapsed_data %>% filter(base >= 31982056 & base <= 32035418)

asterisk_coordinates <- c(31994396, 31994623, 31995781, 31996082, 31996085, 31996093, 31996096, 31996098, 31996450, 31996538, 31996543, 31996552, 31996553, 31996613, 31996617, 31997007, 31997464, 31997605, 32027134, 32027361, 32028519, 32028820, 32028823, 32028831, 32028834, 32028836, 32029188, 32029276, 32029281, 32029290, 32029291, 32029351, 32029355, 32029745, 32030202, 32030343)

# Create the plot using the filtered data
figure <- ggplot(data_filtered, aes(x = base, y = mean_depth)) +
  geom_line(size = 1.25) +  # Apply linewidth 1.25 for all lines
  xlab("Position on Chromosome 6 (Mb)") +  # Update x-axis label to indicate Mb
  ylab("Mean Coverage Depth") + 
  theme(panel.background = element_blank()) +
  theme(panel.border = element_blank(), 
        axis.line = element_line(colour = "black", size = 1), 
        axis.text.x = element_text(size = 16, angle = 45, hjust = 1), 
        axis.text.y = element_text(size = 16), 
        axis.title = element_text(size = 16, face = "bold")) + 
  scale_color_manual(values = sample_colors) +  # Apply the color scale
  scale_x_continuous(
    labels = function(x) sprintf("%.2f", x / 1e6),   # Convert bases to Mb and format with 2 decimal places
    breaks = pretty(data_filtered$base, n = 6)   # Generate 6 tick marks on the x-axis
  ) +
  geom_segment(aes(x = 31982056, y = 299.5, xend = 32002681, yend = 300.5), size = 1, color = "black") +
  annotate("text", x = mean(c(31982056, 32002681)), y = 315, label = "C4A", color = "black", vjust = 0, size = 4) +
  geom_segment(aes(x = 32014794, y = 299.5, xend = 32035418, yend = 300.5), size = 1, color = "black") +
  annotate("text", x = mean(c(32014794, 32035418)), y = 315, label = "C4B", color = "black", vjust = 0, size = 4) + 
  geom_point(data = tibble(x = asterisk_coordinates, y = 330), 
             aes(x = x, y = y), 
             shape = 8, color = "red", size = 3)


# Annotate exons from bed_data
for (i in 1:nrow(bed_data)) {
  figure <- figure + annotate("rect",
                              xmin = bed_data$start[i], 
                              xmax = bed_data$end[i], 
                              ymin = min(290),  
                              ymax = max(310),  
                              fill = "black", color = "black", size = 0.075)
}

figure

ggsave(filename = "c4.pdf", plot = figure, width=169, units = "mm")
